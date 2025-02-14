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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(1.8774348) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(2.0533419) q[1];
sx q[1];
rz(6.2702141) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94840996) q[0];
sx q[0];
rz(-1.4140861) q[0];
sx q[0];
rz(0.47845263) q[0];
x q[1];
rz(1.6799203) q[2];
sx q[2];
rz(-1.732601) q[2];
sx q[2];
rz(-1.2798123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56264191) q[1];
sx q[1];
rz(-1.2516252) q[1];
sx q[1];
rz(-0.36960543) q[1];
x q[2];
rz(0.58297021) q[3];
sx q[3];
rz(-1.7503925) q[3];
sx q[3];
rz(-2.666134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(3.0461779) q[2];
rz(2.6573507) q[3];
sx q[3];
rz(-2.3384194) q[3];
sx q[3];
rz(-1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5937623) q[0];
sx q[0];
rz(-1.7050803) q[0];
sx q[0];
rz(-2.4875212) q[0];
rz(1.7895128) q[1];
sx q[1];
rz(-0.65579954) q[1];
sx q[1];
rz(0.10993122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2034775) q[0];
sx q[0];
rz(-3.118096) q[0];
sx q[0];
rz(2.1142939) q[0];
rz(-pi) q[1];
rz(-0.45808001) q[2];
sx q[2];
rz(-2.1501679) q[2];
sx q[2];
rz(-1.9922076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0866894) q[1];
sx q[1];
rz(-0.48799054) q[1];
sx q[1];
rz(-2.7665374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7450493) q[3];
sx q[3];
rz(-1.3351044) q[3];
sx q[3];
rz(1.1514219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8249417) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(-1.4671154) q[2];
rz(1.0351099) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(-1.3181814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255945) q[0];
sx q[0];
rz(-0.23402973) q[0];
sx q[0];
rz(2.3509534) q[0];
rz(2.3979893) q[1];
sx q[1];
rz(-1.598282) q[1];
sx q[1];
rz(0.57949439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42946896) q[0];
sx q[0];
rz(-0.83659029) q[0];
sx q[0];
rz(0.81487687) q[0];
rz(2.2966301) q[2];
sx q[2];
rz(-0.38533346) q[2];
sx q[2];
rz(-2.8266738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5011141) q[1];
sx q[1];
rz(-1.5334629) q[1];
sx q[1];
rz(1.8289315) q[1];
rz(-pi) q[2];
rz(2.0602588) q[3];
sx q[3];
rz(-2.0474985) q[3];
sx q[3];
rz(-0.34927327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4641331) q[2];
sx q[2];
rz(-2.3049998) q[2];
sx q[2];
rz(0.38132384) q[2];
rz(2.2903806) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(-0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6094991) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(0.36439782) q[0];
rz(2.9413307) q[1];
sx q[1];
rz(-1.4118782) q[1];
sx q[1];
rz(0.099954896) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1276949) q[0];
sx q[0];
rz(-1.5911926) q[0];
sx q[0];
rz(1.4958032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9995297) q[2];
sx q[2];
rz(-1.4754646) q[2];
sx q[2];
rz(-0.39458654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2435088) q[1];
sx q[1];
rz(-1.5931411) q[1];
sx q[1];
rz(-3.0490655) q[1];
rz(0.95619802) q[3];
sx q[3];
rz(-1.0805939) q[3];
sx q[3];
rz(-2.4407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.072307) q[2];
sx q[2];
rz(-1.0682718) q[2];
sx q[2];
rz(-0.11492534) q[2];
rz(-1.140444) q[3];
sx q[3];
rz(-1.8585669) q[3];
sx q[3];
rz(-2.1663402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045778) q[0];
sx q[0];
rz(-0.95252043) q[0];
sx q[0];
rz(0.17247795) q[0];
rz(1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(0.15484658) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.895899) q[0];
sx q[0];
rz(-1.2749199) q[0];
sx q[0];
rz(-2.9135743) q[0];
x q[1];
rz(1.8932976) q[2];
sx q[2];
rz(-2.4716544) q[2];
sx q[2];
rz(2.8084286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78945827) q[1];
sx q[1];
rz(-0.43757868) q[1];
sx q[1];
rz(2.7960013) q[1];
rz(0.17821594) q[3];
sx q[3];
rz(-1.4608188) q[3];
sx q[3];
rz(2.4770155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0871206) q[2];
sx q[2];
rz(-0.27721578) q[2];
sx q[2];
rz(-2.7692914) q[2];
rz(-2.7436658) q[3];
sx q[3];
rz(-1.1436661) q[3];
sx q[3];
rz(0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.953124) q[0];
sx q[0];
rz(-2.5690881) q[0];
sx q[0];
rz(1.5931801) q[0];
rz(-2.7717223) q[1];
sx q[1];
rz(-1.3984171) q[1];
sx q[1];
rz(2.5028548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67154449) q[0];
sx q[0];
rz(-0.28451583) q[0];
sx q[0];
rz(-0.99308641) q[0];
rz(-pi) q[1];
rz(1.2296835) q[2];
sx q[2];
rz(-2.5899867) q[2];
sx q[2];
rz(2.6831339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7914541) q[1];
sx q[1];
rz(-1.119507) q[1];
sx q[1];
rz(-0.43124388) q[1];
rz(-pi) q[2];
rz(-0.5587479) q[3];
sx q[3];
rz(-2.2674446) q[3];
sx q[3];
rz(-3.0697256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0729735) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(0.2529141) q[2];
rz(-0.84826338) q[3];
sx q[3];
rz(-1.4784644) q[3];
sx q[3];
rz(3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0406168) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(0.14895359) q[0];
rz(-1.3111929) q[1];
sx q[1];
rz(-1.7313749) q[1];
sx q[1];
rz(3.0874918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9210756) q[0];
sx q[0];
rz(-1.0722677) q[0];
sx q[0];
rz(-1.4709298) q[0];
rz(-pi) q[1];
rz(0.80246822) q[2];
sx q[2];
rz(-1.6020498) q[2];
sx q[2];
rz(-2.5078338) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91840345) q[1];
sx q[1];
rz(-1.4610054) q[1];
sx q[1];
rz(-2.4415183) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1182958) q[3];
sx q[3];
rz(-1.0956278) q[3];
sx q[3];
rz(3.0123229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80186239) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(2.1745963) q[2];
rz(-2.4407834) q[3];
sx q[3];
rz(-0.98027027) q[3];
sx q[3];
rz(-2.7907659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0130149) q[0];
sx q[0];
rz(-1.9676493) q[0];
sx q[0];
rz(1.5933734) q[0];
rz(2.7633527) q[1];
sx q[1];
rz(-0.75235569) q[1];
sx q[1];
rz(-0.38984782) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8075832) q[0];
sx q[0];
rz(-0.54348677) q[0];
sx q[0];
rz(2.7568211) q[0];
rz(-pi) q[1];
rz(-1.8748449) q[2];
sx q[2];
rz(-1.5549193) q[2];
sx q[2];
rz(2.7823663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9222103) q[1];
sx q[1];
rz(-1.2212906) q[1];
sx q[1];
rz(0.52095682) q[1];
rz(-pi) q[2];
rz(2.8865783) q[3];
sx q[3];
rz(-1.3552291) q[3];
sx q[3];
rz(2.4795585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0924015) q[2];
sx q[2];
rz(-1.1213877) q[2];
sx q[2];
rz(-1.258705) q[2];
rz(0.24426584) q[3];
sx q[3];
rz(-1.786307) q[3];
sx q[3];
rz(-2.145483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8213537) q[0];
sx q[0];
rz(-1.9122253) q[0];
sx q[0];
rz(1.7852596) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(-2.70539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9205775) q[0];
sx q[0];
rz(-1.2599143) q[0];
sx q[0];
rz(-2.6665844) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8955375) q[2];
sx q[2];
rz(-0.51273275) q[2];
sx q[2];
rz(-2.2487179) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0862941) q[1];
sx q[1];
rz(-0.82189593) q[1];
sx q[1];
rz(0.31676745) q[1];
rz(-2.7049191) q[3];
sx q[3];
rz(-1.5092351) q[3];
sx q[3];
rz(-1.7583282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.977481) q[2];
sx q[2];
rz(2.5980921) q[2];
rz(0.049526878) q[3];
sx q[3];
rz(-1.4855569) q[3];
sx q[3];
rz(1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108342) q[0];
sx q[0];
rz(-3.0093091) q[0];
sx q[0];
rz(-0.079205967) q[0];
rz(-0.40009701) q[1];
sx q[1];
rz(-0.85405093) q[1];
sx q[1];
rz(-1.0265464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745747) q[0];
sx q[0];
rz(-1.3368784) q[0];
sx q[0];
rz(-0.83252711) q[0];
x q[1];
rz(-0.36621771) q[2];
sx q[2];
rz(-1.6596268) q[2];
sx q[2];
rz(-0.41761145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0594008) q[1];
sx q[1];
rz(-0.33723661) q[1];
sx q[1];
rz(2.853501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7907294) q[3];
sx q[3];
rz(-1.4961582) q[3];
sx q[3];
rz(0.78001744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93572271) q[2];
sx q[2];
rz(-1.3047855) q[2];
sx q[2];
rz(1.4523466) q[2];
rz(-2.3116889) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5230539) q[0];
sx q[0];
rz(-1.9970311) q[0];
sx q[0];
rz(1.507623) q[0];
rz(1.2670831) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(0.98536678) q[2];
sx q[2];
rz(-1.1078688) q[2];
sx q[2];
rz(-2.3175793) q[2];
rz(-1.9125562) q[3];
sx q[3];
rz(-1.0377586) q[3];
sx q[3];
rz(0.23317045) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
