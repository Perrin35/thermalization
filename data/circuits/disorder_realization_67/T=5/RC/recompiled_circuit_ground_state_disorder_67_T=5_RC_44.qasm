OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23842731) q[0];
sx q[0];
rz(-2.981346) q[0];
sx q[0];
rz(1.2877553) q[0];
rz(2.2924478) q[1];
sx q[1];
rz(-1.949911) q[1];
sx q[1];
rz(-0.75872672) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903126) q[0];
sx q[0];
rz(-2.7669124) q[0];
sx q[0];
rz(2.4640342) q[0];
rz(-2.8992462) q[2];
sx q[2];
rz(-2.0449315) q[2];
sx q[2];
rz(0.15210064) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.030144) q[1];
sx q[1];
rz(-1.2093102) q[1];
sx q[1];
rz(0.8059146) q[1];
rz(-pi) q[2];
rz(-2.4843744) q[3];
sx q[3];
rz(-0.94122488) q[3];
sx q[3];
rz(0.30318015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3005001) q[2];
sx q[2];
rz(-0.79885834) q[2];
sx q[2];
rz(2.63499) q[2];
rz(-1.6182342) q[3];
sx q[3];
rz(-0.62464276) q[3];
sx q[3];
rz(0.055559572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27187207) q[0];
sx q[0];
rz(-0.46590081) q[0];
sx q[0];
rz(0.063657612) q[0];
rz(1.9438538) q[1];
sx q[1];
rz(-1.627219) q[1];
sx q[1];
rz(2.1228085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2301529) q[0];
sx q[0];
rz(-2.8960315) q[0];
sx q[0];
rz(-2.5018238) q[0];
rz(-2.8009861) q[2];
sx q[2];
rz(-1.34769) q[2];
sx q[2];
rz(0.67137137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19551046) q[1];
sx q[1];
rz(-1.2942593) q[1];
sx q[1];
rz(-0.72317381) q[1];
rz(-pi) q[2];
rz(-1.7131931) q[3];
sx q[3];
rz(-0.75800688) q[3];
sx q[3];
rz(0.49331323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7055052) q[2];
sx q[2];
rz(-2.5477396) q[2];
sx q[2];
rz(2.2994821) q[2];
rz(-1.0574794) q[3];
sx q[3];
rz(-2.2130241) q[3];
sx q[3];
rz(-1.9717982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24404003) q[0];
sx q[0];
rz(-2.121448) q[0];
sx q[0];
rz(1.9496339) q[0];
rz(-0.92487088) q[1];
sx q[1];
rz(-0.7083188) q[1];
sx q[1];
rz(-1.8398197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90067139) q[0];
sx q[0];
rz(-2.0427454) q[0];
sx q[0];
rz(-0.94743418) q[0];
x q[1];
rz(1.5315751) q[2];
sx q[2];
rz(-1.7698145) q[2];
sx q[2];
rz(-1.4250371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1731378) q[1];
sx q[1];
rz(-2.0357642) q[1];
sx q[1];
rz(1.370621) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25721217) q[3];
sx q[3];
rz(-1.9055771) q[3];
sx q[3];
rz(-2.0906665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3785582) q[2];
sx q[2];
rz(-0.083746567) q[2];
sx q[2];
rz(0.087510022) q[2];
rz(1.8267501) q[3];
sx q[3];
rz(-2.0851236) q[3];
sx q[3];
rz(0.99353138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4104376) q[0];
sx q[0];
rz(-1.9696099) q[0];
sx q[0];
rz(-2.3325969) q[0];
rz(-1.0357098) q[1];
sx q[1];
rz(-0.46329841) q[1];
sx q[1];
rz(1.0332003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3147723) q[0];
sx q[0];
rz(-2.2345532) q[0];
sx q[0];
rz(-0.04962595) q[0];
rz(-pi) q[1];
rz(-0.43605752) q[2];
sx q[2];
rz(-1.0761217) q[2];
sx q[2];
rz(0.058017284) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82419187) q[1];
sx q[1];
rz(-0.5017952) q[1];
sx q[1];
rz(0.3806033) q[1];
rz(-pi) q[2];
rz(-1.0907034) q[3];
sx q[3];
rz(-1.0786782) q[3];
sx q[3];
rz(1.1034213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0674627) q[2];
sx q[2];
rz(-2.9758657) q[2];
sx q[2];
rz(-1.5376252) q[2];
rz(1.4175203) q[3];
sx q[3];
rz(-1.1105024) q[3];
sx q[3];
rz(1.0874776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23987016) q[0];
sx q[0];
rz(-0.19342315) q[0];
sx q[0];
rz(-1.9523917) q[0];
rz(-2.4173648) q[1];
sx q[1];
rz(-1.291357) q[1];
sx q[1];
rz(1.474818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4714469) q[0];
sx q[0];
rz(-1.0751123) q[0];
sx q[0];
rz(-1.0814654) q[0];
rz(-pi) q[1];
rz(-1.7314265) q[2];
sx q[2];
rz(-2.361627) q[2];
sx q[2];
rz(0.51973625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96383038) q[1];
sx q[1];
rz(-0.46198598) q[1];
sx q[1];
rz(-2.3038008) q[1];
rz(-pi) q[2];
rz(1.0737562) q[3];
sx q[3];
rz(-1.9816363) q[3];
sx q[3];
rz(0.88676329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7597947) q[2];
sx q[2];
rz(-0.90961027) q[2];
sx q[2];
rz(-2.6749715) q[2];
rz(1.7032547) q[3];
sx q[3];
rz(-1.8432901) q[3];
sx q[3];
rz(2.2013825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71317116) q[0];
sx q[0];
rz(-2.3155825) q[0];
sx q[0];
rz(0.010350479) q[0];
rz(-1.6185282) q[1];
sx q[1];
rz(-1.7489988) q[1];
sx q[1];
rz(1.7656322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3165163) q[0];
sx q[0];
rz(-0.42046558) q[0];
sx q[0];
rz(-2.1964873) q[0];
rz(-2.544246) q[2];
sx q[2];
rz(-0.22946363) q[2];
sx q[2];
rz(-0.84537431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.350647) q[1];
sx q[1];
rz(-0.53574359) q[1];
sx q[1];
rz(-0.36046268) q[1];
rz(-pi) q[2];
rz(2.2155432) q[3];
sx q[3];
rz(-0.9315486) q[3];
sx q[3];
rz(3.002142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7213664) q[2];
sx q[2];
rz(-1.0488291) q[2];
sx q[2];
rz(-2.4143207) q[2];
rz(-2.6266802) q[3];
sx q[3];
rz(-1.8102831) q[3];
sx q[3];
rz(-1.7236727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10776831) q[0];
sx q[0];
rz(-1.533778) q[0];
sx q[0];
rz(-0.4069826) q[0];
rz(0.96626967) q[1];
sx q[1];
rz(-0.85710183) q[1];
sx q[1];
rz(3.0028717) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4086558) q[0];
sx q[0];
rz(-1.4545024) q[0];
sx q[0];
rz(3.1058806) q[0];
x q[1];
rz(2.474276) q[2];
sx q[2];
rz(-0.34798589) q[2];
sx q[2];
rz(1.1225357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87200621) q[1];
sx q[1];
rz(-1.9161738) q[1];
sx q[1];
rz(-2.7589382) q[1];
rz(-pi) q[2];
rz(-1.2267031) q[3];
sx q[3];
rz(-1.4315245) q[3];
sx q[3];
rz(-0.95668281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2762642) q[2];
sx q[2];
rz(-1.3778957) q[2];
sx q[2];
rz(2.830937) q[2];
rz(-1.3079414) q[3];
sx q[3];
rz(-1.6303948) q[3];
sx q[3];
rz(2.7697897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0528316) q[0];
sx q[0];
rz(-2.796687) q[0];
sx q[0];
rz(-3.0533691) q[0];
rz(3.131033) q[1];
sx q[1];
rz(-1.5225531) q[1];
sx q[1];
rz(0.47058502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6321736) q[0];
sx q[0];
rz(-2.1996561) q[0];
sx q[0];
rz(0.6412613) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80704524) q[2];
sx q[2];
rz(-1.7084439) q[2];
sx q[2];
rz(-2.0375843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7631665) q[1];
sx q[1];
rz(-1.6554804) q[1];
sx q[1];
rz(1.8316818) q[1];
rz(-pi) q[2];
rz(-3.0017774) q[3];
sx q[3];
rz(-1.2212409) q[3];
sx q[3];
rz(-0.67963723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9528902) q[2];
sx q[2];
rz(-0.95658335) q[2];
sx q[2];
rz(-1.0799705) q[2];
rz(0.88932577) q[3];
sx q[3];
rz(-1.5998799) q[3];
sx q[3];
rz(-1.5573474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66460669) q[0];
sx q[0];
rz(-2.7598858) q[0];
sx q[0];
rz(-0.81740022) q[0];
rz(2.3951702) q[1];
sx q[1];
rz(-2.4669929) q[1];
sx q[1];
rz(0.93961632) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22720756) q[0];
sx q[0];
rz(-3.1244579) q[0];
sx q[0];
rz(1.5987773) q[0];
x q[1];
rz(0.65931658) q[2];
sx q[2];
rz(-2.7020279) q[2];
sx q[2];
rz(-2.2826113) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6691362) q[1];
sx q[1];
rz(-1.8207153) q[1];
sx q[1];
rz(-0.12052287) q[1];
rz(-pi) q[2];
rz(-1.9838748) q[3];
sx q[3];
rz(-1.7677757) q[3];
sx q[3];
rz(0.22991163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32999906) q[2];
sx q[2];
rz(-2.1758695) q[2];
sx q[2];
rz(0.3453556) q[2];
rz(1.5251478) q[3];
sx q[3];
rz(-1.7914146) q[3];
sx q[3];
rz(2.7143872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78373194) q[0];
sx q[0];
rz(-1.1564199) q[0];
sx q[0];
rz(-3.0586437) q[0];
rz(-2.9792765) q[1];
sx q[1];
rz(-1.4444618) q[1];
sx q[1];
rz(0.69581318) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0950553) q[0];
sx q[0];
rz(-2.8588606) q[0];
sx q[0];
rz(0.15034349) q[0];
x q[1];
rz(-1.5478304) q[2];
sx q[2];
rz(-3.0122979) q[2];
sx q[2];
rz(0.33516075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2999794) q[1];
sx q[1];
rz(-2.0149954) q[1];
sx q[1];
rz(1.5919374) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70503321) q[3];
sx q[3];
rz(-1.3516197) q[3];
sx q[3];
rz(2.8104643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7210228) q[2];
sx q[2];
rz(-2.9837954) q[2];
sx q[2];
rz(1.921462) q[2];
rz(-2.5681791) q[3];
sx q[3];
rz(-1.8107332) q[3];
sx q[3];
rz(2.8964608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2777916) q[0];
sx q[0];
rz(-1.5024804) q[0];
sx q[0];
rz(-0.13308751) q[0];
rz(0.0028903891) q[1];
sx q[1];
rz(-0.4160226) q[1];
sx q[1];
rz(-0.70152534) q[1];
rz(2.9691545) q[2];
sx q[2];
rz(-1.3256831) q[2];
sx q[2];
rz(-0.55064461) q[2];
rz(-1.5118128) q[3];
sx q[3];
rz(-2.3702757) q[3];
sx q[3];
rz(2.6487917) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
