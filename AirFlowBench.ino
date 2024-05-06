#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BMP085_U.h>
#include <DHT.h>
#include <LiquidCrystal.h>
#include <math.h>
#define DHTTYPE DHT11
#define DHTPIN_H A1 // DHT sensor placed near to motor head
#define DHTPIN_O A2 // DHT sensor placed near to orifice plate

DHT dhtO(DHTPIN_O, DHTTYPE);
DHT dhtH(DHTPIN_H, DHTTYPE);

LiquidCrystal lcd(8, 9, 4, 5, 6, 7); // Define the pins required by the LCD Keypad Shield

Adafruit_BMP085_Unified bmp = Adafruit_BMP085_Unified(10085);

double DpO_0 = 0;
int a3_0 = 246, i = 0;
bool intake = true;

void setup(void) {
  Serial.begin(115200);
  lcd.begin(16,2); 
  lcd.clear();
  lcd.print(" Air Flow Bench ");
  lcd.setCursor(4,1);    
  lcd.print("starting!");    
  if(!bmp.begin()) { // Initialize the BMP085 sensor (and also the I2C communication)
    Serial.print("Ooops, no BMP085 detected ... Check your wiring or I2C ADDR!");
    lcd.clear();
    lcd.print("    No BMP085");
    lcd.setCursor(4,1);    
    lcd.print("detected!");    
    while(1);
  }
  dhtO.begin(); // DHT11 initialization (SPI communication)
  dhtH.begin(); // DHT11 initialization (SPI communication)
  DpO_0 = readDpSM(200); // Pressure difference across the orifices plate (Pa)
  a3_0 = analogRead(A3); // Pressure difference across the motor head using MPXV5004DP
}

void loop(void) {
  char Orificios[9][10] = {"    1    ", "    2    ", "   1+2   ",
    "   1+3   ", "  1+2+3  ", "  1+3+4  ",
    " 1+2+3+4 ", " 1+2+3+5 ", "1+2+3+4+5"};
  static const double C_v[] = {0.3599326743, 0.6478788137, 1.043804755,
    1.439730697, 2.123602778, 3.023434464,
    3.77929308, 5.03905744, 6.658754474}; // CFM*sqrt[kg/(Pa*m3)] 
  double Q, m, rho, phi, p, DpO, DpH, Cv, pA, pH, rhoH;
  float h = dhtO.readHumidity();
  float T = dhtO.readTemperature();
  phi = h*0.01;

// Get a new BMP085 event
  sensors_event_t event;
  bmp.getEvent(&event);
  if (event.pressure) {
    pA = event.pressure*100.0; // Atmospheric pressure (Pa)
/*
    Serial.print("Pressure:    ");
    Serial.print(event.pressure);
    Serial.println(" hPa");
    float T_BMP; // To store the BMP085 temperature 
    bmp.getTemperature(&T_BMP);
    Serial.print("T_BMP: ");
    Serial.print(T_BMP);
    Serial.println(" °C");
    float seaLevelPressure = SENSORS_PRESSURE_SEALEVELHPA;
    Serial.print("Altitude:    "); 
    Serial.print(bmp.pressureToAltitude(seaLevelPressure,
                                        event.pressure)); 
    Serial.println(" m"); 
*/
  } else {
//    Serial.println("BMP085 sensor error");
    lcd.clear();
    lcd.print(" BMP085 sensor ");
    lcd.setCursor(6,1);    
    lcd.print("error!");    
    return;
  }
  DpH = 4.94*fabs(analogRead(A3)-a3_0); // Pressure difference across the motor head (Pa)
/*  
  Serial.print("DpH: ");
  Serial.print(DpH);
  Serial.println(" Pa");
  Serial.print("DHT Orifice Plate Temperature: ");
  Serial.print(T);
  Serial.println(" °C");
  Serial.print("DHT Orifice Plate Humidity: ");
  Serial.print(h);
  Serial.println(" %");
  Serial.print("DHT Motor Head Temperature: ");
*/
  double T_H = dhtH.readTemperature();
  double h_H = dhtH.readHumidity();
// double phiH = 0.01*h_H;
/*  
  Serial.print(T_H);
  Serial.println(" °C");
  Serial.print("DHT Motor Head Humidity: ");
  Serial.print(h_H);
  Serial.println(" %");
*/
  DpO = fabs(readDpSM(100) - DpO_0); // Pressure difference across the orifices plate (Pa)
/*
  Serial.print("DpO = ");
  Serial.print(DpO);
  Serial.println(" Pa");  
*/
  if (intake) {
    pH = pA - DpH;
    p = pA + DpO;
  } else {
    pH = pA + DpH;
    p  = pA - DpO;
  }
  rho = rhoAr(T,phi,p);
//  rhoH = rhoAr(T_H,phiH,pH);
//  rho = rhoH*p*(T_H+273.15)/(pH*(T+273.15));
/*
  Serial.print("rho = ");
  Serial.print(rho);
  Serial.println(" kg/m3");
*/
  Q = C_v[i]*sqrt(2*DpO/rho); // Air flow (CFM)
  m = rho*0.47194745*Q; // Mass air flow (g/s)
  Serial.println(m, 3);
/*
  Serial.print("m = ");
  Serial.print(m);
  Serial.println(" g/s");
  Serial.println("");
*/    
  delay(2000);
  int LCDKpad = analogRead(A0); // Signal coming from the keypad
  if (LCDKpad < 60) { // Second button from right to left pressed (the first is for reset)
    lcd.setCursor(0,0);  
    lcd.print("rho = ");
    lcd.print(rho,2);
    lcd.print("      ");
    lcd.setCursor(11,0);
    lcd.print("kg/m3");
    lcd.setCursor(0,1);
    lcd.print("pA = ");
    lcd.print(pA,1);
    lcd.print("      ");
    lcd.setCursor(13,1);
    lcd.print(" Pa");
  } else if (LCDKpad < 200) { // Upper button pressed
    if (intake) {
      intake = false;
    } else {
      i++;
      if (i > 8) i = 0;
    }
    Q = C_v[i]*sqrt(2*DpO/rho);
    lcd.setCursor(0,0);  
    lcd.print("Q = ");
    lcd.print(Q,1);
    lcd.print("        ");
    lcd.setCursor(13,0);
    lcd.print("CFM");
    lcd.setCursor(0,1);
    lcd.print("Exhst. ");
    lcd.setCursor(7,1);
    lcd.print(Orificios[i]);
  } else if (LCDKpad < 400) { // Lower button pressed
    if (intake) {
      i++;
      if (i > 8) i = 0;
    } else {
      intake = true;
    }
    Q = C_v[i]*sqrt(2*DpO/rho);
    lcd.setCursor(0,0);  
    lcd.print("Q = ");
    lcd.print(Q,1);
    lcd.print("        ");
    lcd.setCursor(13,0);
    lcd.print("CFM");
    lcd.setCursor(0,1);
    lcd.print("Intake ");
    lcd.setCursor(7,1);
    lcd.print(Orificios[i]);
  } else if (LCDKpad < 600) { // Second button from left to right pressed
    lcd.setCursor(0,0);  
    lcd.print("T = ");
    lcd.print(T,1);
    lcd.print("        ");
    lcd.setCursor(13,0);
    lcd.print(" C ");
    lcd.setCursor(0,1);
    lcd.print("UR = ");
    lcd.print(h,1);
    lcd.print("        ");
    lcd.setCursor(13,1);
    lcd.print(" % ");
  } else if (LCDKpad < 800) { // First button from left to right pressed
    lcd.setCursor(0,0);  
    lcd.print("DpH= ");
    lcd.print(DpH/(9.81*25.4),2);
    lcd.print("     ");
    lcd.setCursor(11,0);
    lcd.print("inH2O");
    lcd.setCursor(0,1);
    lcd.print("T_H= ");
    lcd.print(T_H,1);
    lcd.print("        ");
    lcd.setCursor(14,1);
    lcd.print(" C");
  } else { // No buttons pressed
    lcd.setCursor(0,0);  
    lcd.print("Q = ");
    lcd.print(Q,1);
    lcd.print("        ");
    lcd.setCursor(13,0);
    lcd.print("CFM");
    lcd.setCursor(0,1);
    if (intake) {
      lcd.print("Intake ");
    } else {
    lcd.print("Exhst. ");
    }    
    lcd.setCursor(7,1);
    lcd.print(Orificios[i]);    
  }
}

// Air density calculation as a function of temperatura, relative humidity and pressure

float rhoAr(double T, double phi, double p) {
  double p_ar, p_H2O, p_sat, aux;
  static const double M_ar=0.028964,M_H2O=0.018016,R=8.3145;
  static const double e_s0=610.78, c[]={0.99999683,-0.90826951e-2,0.78736169e-4,
    -0.61117958e-6,0.43884187e-8,-0.29883885e-10,0.21874425e-12,
    -0.17892321e-14,0.11112018e-16,-0.30994571e-19};
//Calculo da pressao de vapor saturado
  aux = c[9];
  for (int i=8; i>=0 ;i--)
    aux = c[i] + T*aux;
  p_sat = e_s0*pow(aux,-8);
//Calculo da pressao parcial de vapor
  p_H2O = phi*p_sat;
//Calculo da pressao parcial de ar seco
  p_ar = p - p_H2O;
//Calculo da massa especifica do ar umido
  return((p_ar*M_ar + p_H2O*M_H2O)/(R*(T+273.15)));
}

// Pressure difference reading from SM7331-BCE-S-500.00-351

#define SM_DP_MAX 500.0
#define SM_DP_MIN -500.0
#define SM_DFS 65535.0
#define SM_ADD 0x6C
#define SM_DP_ADD 0x30

double readDpSM(int nSamples) {
  int16_t dout;
  uint8_t lo_byte;
  int8_t hi_byte;
  double dp = 0.0;

  Wire.beginTransmission(byte(SM_ADD));
  Wire.write(byte(SM_DP_ADD));
  Wire.endTransmission();
  delay(50);
  for (int i=0; i<nSamples; i++) {
      Wire.requestFrom(byte(SM_ADD), 2);
      lo_byte = Wire.read(); // receive low byte as lower 8 bits (unsigned)
      hi_byte = Wire.read(); // receive high byte (signed)
      dout = (hi_byte << 8) | lo_byte; // shift high byte to be high 8 bits and completes the reading with low byte
      dp += (SM_DP_MAX - SM_DP_MIN)/nSamples*double(dout)/SM_DFS;
  }
  return(dp + 0.5*(SM_DP_MAX + SM_DP_MIN)); // pressure difference in Pa (averaged over nSamples samples)
}