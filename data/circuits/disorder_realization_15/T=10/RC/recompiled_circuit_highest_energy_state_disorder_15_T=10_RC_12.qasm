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
rz(1.8279561) q[0];
sx q[0];
rz(-0.067129048) q[0];
sx q[0];
rz(0.022775291) q[0];
rz(1.0442806) q[1];
sx q[1];
rz(-2.2106946) q[1];
sx q[1];
rz(0.013962362) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9305796) q[0];
sx q[0];
rz(-1.530181) q[0];
sx q[0];
rz(0.75675772) q[0];
x q[1];
rz(2.3910782) q[2];
sx q[2];
rz(-2.566411) q[2];
sx q[2];
rz(-1.8668979) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.57052862) q[1];
sx q[1];
rz(-1.7693874) q[1];
sx q[1];
rz(1.3873439) q[1];
rz(-0.45043378) q[3];
sx q[3];
rz(-1.2518468) q[3];
sx q[3];
rz(1.3802841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17249168) q[2];
sx q[2];
rz(-1.3250985) q[2];
sx q[2];
rz(-1.3410404) q[2];
rz(1.4260882) q[3];
sx q[3];
rz(-1.1172349) q[3];
sx q[3];
rz(0.30556998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7476927) q[0];
sx q[0];
rz(-1.1423528) q[0];
sx q[0];
rz(0.72545141) q[0];
rz(-2.2254288) q[1];
sx q[1];
rz(-0.80892816) q[1];
sx q[1];
rz(0.61942548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5529157) q[0];
sx q[0];
rz(-1.1631794) q[0];
sx q[0];
rz(2.1216395) q[0];
rz(-2.393707) q[2];
sx q[2];
rz(-1.1074142) q[2];
sx q[2];
rz(3.0836251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.404285) q[1];
sx q[1];
rz(-2.5661307) q[1];
sx q[1];
rz(2.2810443) q[1];
rz(-pi) q[2];
rz(1.0221865) q[3];
sx q[3];
rz(-2.6319176) q[3];
sx q[3];
rz(1.0773848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3118423) q[2];
sx q[2];
rz(-1.0210911) q[2];
sx q[2];
rz(-0.96168438) q[2];
rz(1.1059777) q[3];
sx q[3];
rz(-2.1092236) q[3];
sx q[3];
rz(2.7963514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85213929) q[0];
sx q[0];
rz(-2.0134605) q[0];
sx q[0];
rz(1.3683251) q[0];
rz(0.013956919) q[1];
sx q[1];
rz(-1.114926) q[1];
sx q[1];
rz(1.7272635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089979261) q[0];
sx q[0];
rz(-2.8169301) q[0];
sx q[0];
rz(1.2847631) q[0];
x q[1];
rz(1.3886202) q[2];
sx q[2];
rz(-2.7668556) q[2];
sx q[2];
rz(-2.665525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77905475) q[1];
sx q[1];
rz(-2.0749054) q[1];
sx q[1];
rz(-2.6466531) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8295123) q[3];
sx q[3];
rz(-1.95041) q[3];
sx q[3];
rz(-0.9791475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3783375) q[2];
sx q[2];
rz(-1.9941284) q[2];
sx q[2];
rz(-2.9249127) q[2];
rz(2.2583151) q[3];
sx q[3];
rz(-2.7933385) q[3];
sx q[3];
rz(-2.0047552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9388409) q[0];
sx q[0];
rz(-1.0067679) q[0];
sx q[0];
rz(2.3034565) q[0];
rz(2.5996767) q[1];
sx q[1];
rz(-0.37626615) q[1];
sx q[1];
rz(-3.0964877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7411557) q[0];
sx q[0];
rz(-0.61624762) q[0];
sx q[0];
rz(0.23016302) q[0];
rz(-pi) q[1];
rz(1.0884645) q[2];
sx q[2];
rz(-1.8016029) q[2];
sx q[2];
rz(2.7969691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.038626898) q[1];
sx q[1];
rz(-0.47768394) q[1];
sx q[1];
rz(-1.4250907) q[1];
rz(-pi) q[2];
rz(-0.18121775) q[3];
sx q[3];
rz(-0.50091302) q[3];
sx q[3];
rz(1.25327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0185467) q[2];
sx q[2];
rz(-0.16877731) q[2];
sx q[2];
rz(0.65400845) q[2];
rz(2.5096014) q[3];
sx q[3];
rz(-1.4258823) q[3];
sx q[3];
rz(0.33647195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80394799) q[0];
sx q[0];
rz(-2.4027282) q[0];
sx q[0];
rz(-1.3496572) q[0];
rz(-2.3140287) q[1];
sx q[1];
rz(-0.63500985) q[1];
sx q[1];
rz(0.48447022) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.890936) q[0];
sx q[0];
rz(-1.93205) q[0];
sx q[0];
rz(2.7879232) q[0];
rz(-pi) q[1];
rz(-1.7155649) q[2];
sx q[2];
rz(-2.5826404) q[2];
sx q[2];
rz(-1.1048855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4836241) q[1];
sx q[1];
rz(-1.6977024) q[1];
sx q[1];
rz(0.42327858) q[1];
x q[2];
rz(0.87501672) q[3];
sx q[3];
rz(-1.2051123) q[3];
sx q[3];
rz(0.48236267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1822002) q[2];
sx q[2];
rz(-2.0727938) q[2];
sx q[2];
rz(0.63641316) q[2];
rz(-0.098879769) q[3];
sx q[3];
rz(-1.5560047) q[3];
sx q[3];
rz(-1.4720565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32475489) q[0];
sx q[0];
rz(-0.88961283) q[0];
sx q[0];
rz(1.1725934) q[0];
rz(1.532754) q[1];
sx q[1];
rz(-2.1295261) q[1];
sx q[1];
rz(-3.1408659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9726192) q[0];
sx q[0];
rz(-1.0797636) q[0];
sx q[0];
rz(0.48788496) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9209599) q[2];
sx q[2];
rz(-1.4743525) q[2];
sx q[2];
rz(-0.19410832) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.128617) q[1];
sx q[1];
rz(-0.3518577) q[1];
sx q[1];
rz(0.33280876) q[1];
rz(-0.29720593) q[3];
sx q[3];
rz(-2.5416221) q[3];
sx q[3];
rz(1.456049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26649228) q[2];
sx q[2];
rz(-2.2537587) q[2];
sx q[2];
rz(-0.20981851) q[2];
rz(2.3328414) q[3];
sx q[3];
rz(-0.54976141) q[3];
sx q[3];
rz(-0.89890283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6784191) q[0];
sx q[0];
rz(-1.9371978) q[0];
sx q[0];
rz(0.11209442) q[0];
rz(3.0922999) q[1];
sx q[1];
rz(-2.6105328) q[1];
sx q[1];
rz(1.3370399) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5264915) q[0];
sx q[0];
rz(-1.3630629) q[0];
sx q[0];
rz(-2.2206942) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0626917) q[2];
sx q[2];
rz(-1.4586978) q[2];
sx q[2];
rz(-2.1014121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2087792) q[1];
sx q[1];
rz(-1.5808388) q[1];
sx q[1];
rz(0.34681635) q[1];
x q[2];
rz(-1.2483761) q[3];
sx q[3];
rz(-0.60203505) q[3];
sx q[3];
rz(2.4436434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8423395) q[2];
sx q[2];
rz(-2.3167593) q[2];
sx q[2];
rz(0.71005589) q[2];
rz(3.1402785) q[3];
sx q[3];
rz(-2.240286) q[3];
sx q[3];
rz(0.60091758) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173219) q[0];
sx q[0];
rz(-1.3794427) q[0];
sx q[0];
rz(-0.58260179) q[0];
rz(-2.461589) q[1];
sx q[1];
rz(-0.99907196) q[1];
sx q[1];
rz(-1.2635788) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.82915) q[0];
sx q[0];
rz(-2.0264605) q[0];
sx q[0];
rz(1.8862886) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81496691) q[2];
sx q[2];
rz(-2.8475347) q[2];
sx q[2];
rz(-0.61758274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78681011) q[1];
sx q[1];
rz(-0.92325961) q[1];
sx q[1];
rz(1.3446273) q[1];
rz(0.42383343) q[3];
sx q[3];
rz(-1.5409711) q[3];
sx q[3];
rz(2.3526997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6363643) q[2];
sx q[2];
rz(-1.4437081) q[2];
sx q[2];
rz(3.0155724) q[2];
rz(1.5382918) q[3];
sx q[3];
rz(-0.77682972) q[3];
sx q[3];
rz(-2.7974424) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4890323) q[0];
sx q[0];
rz(-1.6283988) q[0];
sx q[0];
rz(-1.9788347) q[0];
rz(-1.8310422) q[1];
sx q[1];
rz(-2.4030011) q[1];
sx q[1];
rz(1.383925) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726295) q[0];
sx q[0];
rz(-1.9101189) q[0];
sx q[0];
rz(3.0940381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53015401) q[2];
sx q[2];
rz(-2.9994476) q[2];
sx q[2];
rz(1.7696385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0392929) q[1];
sx q[1];
rz(-2.1792534) q[1];
sx q[1];
rz(3.0000172) q[1];
rz(-pi) q[2];
x q[2];
rz(0.042785809) q[3];
sx q[3];
rz(-1.7093911) q[3];
sx q[3];
rz(2.3530172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80013529) q[2];
sx q[2];
rz(-1.8693482) q[2];
sx q[2];
rz(-2.6696491) q[2];
rz(-0.74884549) q[3];
sx q[3];
rz(-2.2180836) q[3];
sx q[3];
rz(-0.81764618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5188468) q[0];
sx q[0];
rz(-2.2846344) q[0];
sx q[0];
rz(-2.6527606) q[0];
rz(-2.7986616) q[1];
sx q[1];
rz(-2.2551426) q[1];
sx q[1];
rz(1.0541281) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3801479) q[0];
sx q[0];
rz(-1.6223546) q[0];
sx q[0];
rz(0.037550302) q[0];
rz(-pi) q[1];
rz(2.6661886) q[2];
sx q[2];
rz(-2.4717071) q[2];
sx q[2];
rz(2.6500882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78275241) q[1];
sx q[1];
rz(-0.61780518) q[1];
sx q[1];
rz(-1.0516404) q[1];
rz(-3.0929186) q[3];
sx q[3];
rz(-2.0423539) q[3];
sx q[3];
rz(-0.31192985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1663345) q[2];
sx q[2];
rz(-2.9322093) q[2];
sx q[2];
rz(-3.0496822) q[2];
rz(-2.6722243) q[3];
sx q[3];
rz(-2.2769603) q[3];
sx q[3];
rz(0.11274591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54409201) q[0];
sx q[0];
rz(-1.8202029) q[0];
sx q[0];
rz(1.2910917) q[0];
rz(-2.3757833) q[1];
sx q[1];
rz(-1.3594834) q[1];
sx q[1];
rz(1.8654738) q[1];
rz(-1.932939) q[2];
sx q[2];
rz(-1.3315143) q[2];
sx q[2];
rz(-1.3523921) q[2];
rz(1.3113931) q[3];
sx q[3];
rz(-1.048199) q[3];
sx q[3];
rz(-0.38182624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
