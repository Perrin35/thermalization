OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3408605) q[0];
sx q[0];
rz(-2.1600318) q[0];
sx q[0];
rz(-3.0032934) q[0];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.9041833) q[2];
sx q[2];
rz(-0.29793973) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.032174234) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(-1.4515619) q[1];
rz(1.4708038) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-0.28764763) q[2];
rz(-1.7488165) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47050414) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(-1.6862306) q[0];
rz(-pi) q[1];
rz(2.4096476) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(1.0674455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8375081) q[1];
sx q[1];
rz(-1.6581074) q[1];
sx q[1];
rz(-2.0259894) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3489477) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-1.0822901) q[2];
rz(1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2079175) q[0];
sx q[0];
rz(-1.9662004) q[0];
sx q[0];
rz(1.0008706) q[0];
rz(-pi) q[1];
rz(0.47282131) q[2];
sx q[2];
rz(-0.46337767) q[2];
sx q[2];
rz(-0.18714999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5971165) q[1];
sx q[1];
rz(-0.58864486) q[1];
sx q[1];
rz(-2.4878923) q[1];
rz(-pi) q[2];
rz(0.28762443) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(2.6401816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(-0.71980643) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(0.71587193) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-2.148927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(0.53702766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56447345) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(1.0587453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.042295481) q[1];
sx q[1];
rz(-0.18810939) q[1];
sx q[1];
rz(0.52109615) q[1];
x q[2];
rz(0.4977787) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(-1.9073515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0585534) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(1.7129718) q[2];
rz(-0.84609091) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(2.3748421) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(-2.7899182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(2.48824) q[0];
rz(-pi) q[1];
rz(-2.211116) q[2];
sx q[2];
rz(-2.3978007) q[2];
sx q[2];
rz(3.004068) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.60492086) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(-0.037996304) q[1];
rz(-pi) q[2];
rz(-2.8204927) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(2.0146446) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-2.0419962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27605155) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(-1.0827351) q[0];
x q[1];
rz(-0.18892388) q[2];
sx q[2];
rz(-2.0940229) q[2];
sx q[2];
rz(-0.97666937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3179143) q[1];
sx q[1];
rz(-2.4412529) q[1];
sx q[1];
rz(-2.8631696) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33575616) q[3];
sx q[3];
rz(-0.81695518) q[3];
sx q[3];
rz(-2.9591065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(0.28298322) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.5420087) q[3];
sx q[3];
rz(2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(-0.15730102) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(1.69676) q[0];
rz(2.3072725) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(-1.8404567) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7423319) q[1];
sx q[1];
rz(-2.4943588) q[1];
sx q[1];
rz(-1.8982235) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8623452) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1823696) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-2.9390745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-0.29327926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54942185) q[0];
sx q[0];
rz(-0.54245078) q[0];
sx q[0];
rz(-0.011499238) q[0];
rz(-pi) q[1];
rz(-1.3515527) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(3.0041681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9353232) q[1];
sx q[1];
rz(-2.2895797) q[1];
sx q[1];
rz(-2.1410336) q[1];
rz(-pi) q[2];
rz(-1.1612438) q[3];
sx q[3];
rz(-2.000994) q[3];
sx q[3];
rz(-2.6524515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(0.088767178) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(1.3508505) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1254743) q[0];
sx q[0];
rz(-3.0258412) q[0];
sx q[0];
rz(0.8158169) q[0];
rz(-pi) q[1];
rz(0.10599372) q[2];
sx q[2];
rz(-2.4786737) q[2];
sx q[2];
rz(1.7352599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6752154) q[1];
sx q[1];
rz(-1.5545168) q[1];
sx q[1];
rz(-1.1364163) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1397347) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(-1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(-0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7923332) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(2.7381998) q[0];
rz(3.0868297) q[2];
sx q[2];
rz(-0.4224531) q[2];
sx q[2];
rz(1.7516608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.032220275) q[1];
sx q[1];
rz(-1.6175555) q[1];
sx q[1];
rz(0.95581518) q[1];
rz(-pi) q[2];
rz(-2.9986936) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(-2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(2.4492241) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(2.6295173) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(-0.84531534) q[2];
sx q[2];
rz(-0.39388638) q[2];
sx q[2];
rz(-2.3245927) q[2];
rz(-0.11733304) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
