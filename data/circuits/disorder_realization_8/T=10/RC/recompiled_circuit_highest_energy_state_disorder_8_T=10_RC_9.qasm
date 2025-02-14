OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(-3.0515766) q[1];
sx q[1];
rz(-0.4740735) q[1];
sx q[1];
rz(0.83516821) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1195364) q[0];
sx q[0];
rz(-1.8827264) q[0];
sx q[0];
rz(-2.3989912) q[0];
x q[1];
rz(1.2537755) q[2];
sx q[2];
rz(-0.95404139) q[2];
sx q[2];
rz(0.9483017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0402643) q[1];
sx q[1];
rz(-0.44875408) q[1];
sx q[1];
rz(-0.59228102) q[1];
rz(-pi) q[2];
rz(0.6938758) q[3];
sx q[3];
rz(-0.78842794) q[3];
sx q[3];
rz(0.83886787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-2.2455402) q[2];
sx q[2];
rz(-0.33207616) q[2];
rz(2.8927228) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(-2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2700972) q[0];
sx q[0];
rz(-2.8506554) q[0];
sx q[0];
rz(2.6089597) q[0];
rz(-0.10781413) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(-0.32049387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4491534) q[0];
sx q[0];
rz(-2.4863613) q[0];
sx q[0];
rz(-1.6754608) q[0];
rz(1.8024496) q[2];
sx q[2];
rz(-1.7391053) q[2];
sx q[2];
rz(1.2141808) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1793666) q[1];
sx q[1];
rz(-1.4186064) q[1];
sx q[1];
rz(-0.16140143) q[1];
rz(-2.8122607) q[3];
sx q[3];
rz(-1.3668225) q[3];
sx q[3];
rz(-2.4402809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8104441) q[2];
sx q[2];
rz(-0.30511567) q[2];
sx q[2];
rz(1.6395052) q[2];
rz(2.8034927) q[3];
sx q[3];
rz(-2.2564042) q[3];
sx q[3];
rz(-0.52687183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6239887) q[0];
sx q[0];
rz(-1.9094587) q[0];
sx q[0];
rz(0.35743085) q[0];
rz(2.7492211) q[1];
sx q[1];
rz(-2.3452499) q[1];
sx q[1];
rz(1.2145112) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314963) q[0];
sx q[0];
rz(-1.4292681) q[0];
sx q[0];
rz(-0.017472762) q[0];
rz(-pi) q[1];
rz(0.54698555) q[2];
sx q[2];
rz(-1.437709) q[2];
sx q[2];
rz(-2.2625429) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9704772) q[1];
sx q[1];
rz(-1.5678645) q[1];
sx q[1];
rz(1.9582002) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1585095) q[3];
sx q[3];
rz(-1.6020892) q[3];
sx q[3];
rz(-0.33558095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3452722) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(-0.39682445) q[2];
rz(1.2567629) q[3];
sx q[3];
rz(-1.4182914) q[3];
sx q[3];
rz(-2.5313012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20763718) q[0];
sx q[0];
rz(-1.7455245) q[0];
sx q[0];
rz(1.9130094) q[0];
rz(-1.2127016) q[1];
sx q[1];
rz(-0.96114254) q[1];
sx q[1];
rz(-2.520715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69260864) q[0];
sx q[0];
rz(-1.0978069) q[0];
sx q[0];
rz(-0.18758298) q[0];
rz(0.75350301) q[2];
sx q[2];
rz(-2.916159) q[2];
sx q[2];
rz(-1.6627251) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.828307) q[1];
sx q[1];
rz(-1.3125712) q[1];
sx q[1];
rz(-0.26744803) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38264783) q[3];
sx q[3];
rz(-1.4605902) q[3];
sx q[3];
rz(0.73303849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.88317251) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(-0.15667285) q[2];
rz(2.2251718) q[3];
sx q[3];
rz(-1.0333002) q[3];
sx q[3];
rz(-2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857392) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(2.7440985) q[0];
rz(0.022857895) q[1];
sx q[1];
rz(-2.6409179) q[1];
sx q[1];
rz(0.674725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8431339) q[0];
sx q[0];
rz(-1.805691) q[0];
sx q[0];
rz(2.9265334) q[0];
x q[1];
rz(0.12395383) q[2];
sx q[2];
rz(-1.6646241) q[2];
sx q[2];
rz(-1.7403062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9123332) q[1];
sx q[1];
rz(-2.9119171) q[1];
sx q[1];
rz(2.4116184) q[1];
rz(-1.2615023) q[3];
sx q[3];
rz(-1.8320683) q[3];
sx q[3];
rz(-2.3051534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(0.69925365) q[2];
rz(1.3462542) q[3];
sx q[3];
rz(-0.56763879) q[3];
sx q[3];
rz(-2.4301372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72494495) q[0];
sx q[0];
rz(-2.7236433) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(-0.97459546) q[1];
sx q[1];
rz(-1.83788) q[1];
sx q[1];
rz(-1.4385361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33773003) q[0];
sx q[0];
rz(-0.32050214) q[0];
sx q[0];
rz(-0.33945531) q[0];
rz(-pi) q[1];
rz(0.010377093) q[2];
sx q[2];
rz(-1.8266018) q[2];
sx q[2];
rz(2.1880045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4230835) q[1];
sx q[1];
rz(-0.72186493) q[1];
sx q[1];
rz(1.9265429) q[1];
rz(1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(-2.0723267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4633816) q[2];
sx q[2];
rz(-1.7794098) q[2];
sx q[2];
rz(-1.8360651) q[2];
rz(1.3156923) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(-0.8684043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0618133) q[0];
sx q[0];
rz(-2.7257305) q[0];
sx q[0];
rz(1.4965936) q[0];
rz(0.98622259) q[1];
sx q[1];
rz(-1.5602427) q[1];
sx q[1];
rz(2.644002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3517036) q[0];
sx q[0];
rz(-1.6663807) q[0];
sx q[0];
rz(1.315675) q[0];
x q[1];
rz(-0.87734434) q[2];
sx q[2];
rz(-2.3893271) q[2];
sx q[2];
rz(1.8107506) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6141245) q[1];
sx q[1];
rz(-1.9640199) q[1];
sx q[1];
rz(0.16904633) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4769745) q[3];
sx q[3];
rz(-1.854343) q[3];
sx q[3];
rz(-0.16167262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8387973) q[2];
sx q[2];
rz(-1.5747728) q[2];
sx q[2];
rz(2.8509169) q[2];
rz(-2.931328) q[3];
sx q[3];
rz(-0.99431521) q[3];
sx q[3];
rz(2.1073585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5656089) q[0];
sx q[0];
rz(-1.9714332) q[0];
sx q[0];
rz(-1.959311) q[0];
rz(0.15444175) q[1];
sx q[1];
rz(-1.7223822) q[1];
sx q[1];
rz(-0.90528893) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2574999) q[0];
sx q[0];
rz(-1.582889) q[0];
sx q[0];
rz(-1.2980677) q[0];
rz(-pi) q[1];
rz(0.45640517) q[2];
sx q[2];
rz(-1.3564566) q[2];
sx q[2];
rz(-0.53949196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33069995) q[1];
sx q[1];
rz(-2.3726683) q[1];
sx q[1];
rz(0.21067174) q[1];
rz(-0.21568294) q[3];
sx q[3];
rz(-1.9043916) q[3];
sx q[3];
rz(1.622596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0514544) q[2];
sx q[2];
rz(-1.5807512) q[2];
sx q[2];
rz(2.1913989) q[2];
rz(-0.072362445) q[3];
sx q[3];
rz(-1.7642998) q[3];
sx q[3];
rz(-0.70934057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41910928) q[0];
sx q[0];
rz(-2.1868732) q[0];
sx q[0];
rz(-1.0954274) q[0];
rz(-1.0911881) q[1];
sx q[1];
rz(-2.2406816) q[1];
sx q[1];
rz(-0.99517623) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93928775) q[0];
sx q[0];
rz(-1.5065333) q[0];
sx q[0];
rz(0.66575428) q[0];
rz(0.27492304) q[2];
sx q[2];
rz(-1.2464393) q[2];
sx q[2];
rz(2.7664326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6217664) q[1];
sx q[1];
rz(-1.4859746) q[1];
sx q[1];
rz(1.6165401) q[1];
rz(1.770129) q[3];
sx q[3];
rz(-1.3889894) q[3];
sx q[3];
rz(2.0169286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5114078) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(2.4616145) q[2];
rz(1.4558815) q[3];
sx q[3];
rz(-2.2472491) q[3];
sx q[3];
rz(-1.1792012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42530123) q[0];
sx q[0];
rz(-2.20708) q[0];
sx q[0];
rz(2.5262078) q[0];
rz(1.8062493) q[1];
sx q[1];
rz(-1.6067303) q[1];
sx q[1];
rz(-1.3478442) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6123209) q[0];
sx q[0];
rz(-2.7833496) q[0];
sx q[0];
rz(3.0105581) q[0];
x q[1];
rz(2.2632019) q[2];
sx q[2];
rz(-0.88365245) q[2];
sx q[2];
rz(0.85604446) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8948095) q[1];
sx q[1];
rz(-1.9468773) q[1];
sx q[1];
rz(-1.7433107) q[1];
rz(-pi) q[2];
rz(1.4487292) q[3];
sx q[3];
rz(-0.81840912) q[3];
sx q[3];
rz(-0.19644745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24796692) q[2];
sx q[2];
rz(-1.8411571) q[2];
sx q[2];
rz(2.628053) q[2];
rz(2.9945471) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(-1.8769544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2697987) q[0];
sx q[0];
rz(-2.2525621) q[0];
sx q[0];
rz(2.9004108) q[0];
rz(-1.8409894) q[1];
sx q[1];
rz(-2.0722957) q[1];
sx q[1];
rz(3.121079) q[1];
rz(1.9199964) q[2];
sx q[2];
rz(-0.39633718) q[2];
sx q[2];
rz(0.63971165) q[2];
rz(-0.5823395) q[3];
sx q[3];
rz(-2.4262541) q[3];
sx q[3];
rz(-2.2179009) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
