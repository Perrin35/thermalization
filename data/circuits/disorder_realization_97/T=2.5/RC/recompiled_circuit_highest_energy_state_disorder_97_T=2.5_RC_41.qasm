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
rz(-0.66145575) q[0];
sx q[0];
rz(2.8602726) q[0];
sx q[0];
rz(8.2814132) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(-1.3209359) q[1];
sx q[1];
rz(0.015425711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93958873) q[0];
sx q[0];
rz(-2.1103854) q[0];
sx q[0];
rz(-0.57641502) q[0];
rz(-0.77818971) q[2];
sx q[2];
rz(-0.74290771) q[2];
sx q[2];
rz(-1.5404494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.4651067) q[1];
sx q[1];
rz(-2.3576405) q[1];
sx q[1];
rz(-0.98760651) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32415819) q[3];
sx q[3];
rz(-2.5796842) q[3];
sx q[3];
rz(-3.0747385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.780484) q[2];
sx q[2];
rz(-0.0094272308) q[2];
sx q[2];
rz(-1.1040322) q[2];
rz(-2.5822254) q[3];
sx q[3];
rz(-0.95061022) q[3];
sx q[3];
rz(2.2573788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169432) q[0];
sx q[0];
rz(-0.58498061) q[0];
sx q[0];
rz(0.90644932) q[0];
rz(-2.8001884) q[1];
sx q[1];
rz(-0.87541348) q[1];
sx q[1];
rz(0.34814775) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43480647) q[0];
sx q[0];
rz(-0.40074391) q[0];
sx q[0];
rz(-0.67760076) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0013214) q[2];
sx q[2];
rz(-3.0932326) q[2];
sx q[2];
rz(2.0404301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2359008) q[1];
sx q[1];
rz(-2.5254619) q[1];
sx q[1];
rz(1.4409164) q[1];
rz(-2.9112629) q[3];
sx q[3];
rz(-0.84127142) q[3];
sx q[3];
rz(-1.426633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97588313) q[2];
sx q[2];
rz(-2.7996863) q[2];
sx q[2];
rz(-0.085414097) q[2];
rz(0.092770569) q[3];
sx q[3];
rz(-0.86920357) q[3];
sx q[3];
rz(0.87576491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669325) q[0];
sx q[0];
rz(-2.7800738) q[0];
sx q[0];
rz(-0.34459484) q[0];
rz(-1.5498281) q[1];
sx q[1];
rz(-1.4180309) q[1];
sx q[1];
rz(1.6835015) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8710468) q[0];
sx q[0];
rz(-2.4231909) q[0];
sx q[0];
rz(-3.0860391) q[0];
rz(-pi) q[1];
rz(-3.0821716) q[2];
sx q[2];
rz(-2.0869617) q[2];
sx q[2];
rz(-2.8490861) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96839206) q[1];
sx q[1];
rz(-1.3276982) q[1];
sx q[1];
rz(2.9673884) q[1];
x q[2];
rz(2.273748) q[3];
sx q[3];
rz(-0.72753564) q[3];
sx q[3];
rz(-0.65505469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.496326) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(-2.0110896) q[2];
rz(-0.94275236) q[3];
sx q[3];
rz(-1.0470942) q[3];
sx q[3];
rz(-1.9345136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4261674) q[0];
sx q[0];
rz(-1.796145) q[0];
sx q[0];
rz(-1.4581534) q[0];
rz(1.0484877) q[1];
sx q[1];
rz(-1.6494992) q[1];
sx q[1];
rz(-2.6715211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94435197) q[0];
sx q[0];
rz(-1.2908123) q[0];
sx q[0];
rz(-2.6231595) q[0];
rz(-2.0087035) q[2];
sx q[2];
rz(-0.57269579) q[2];
sx q[2];
rz(1.1066135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64420358) q[1];
sx q[1];
rz(-2.4371689) q[1];
sx q[1];
rz(-2.8423487) q[1];
x q[2];
rz(0.87155452) q[3];
sx q[3];
rz(-2.1208753) q[3];
sx q[3];
rz(-3.0015903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5992392) q[2];
sx q[2];
rz(-0.75672823) q[2];
sx q[2];
rz(2.2139464) q[2];
rz(1.2841691) q[3];
sx q[3];
rz(-1.410306) q[3];
sx q[3];
rz(-0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1272142) q[0];
sx q[0];
rz(-2.7365186) q[0];
sx q[0];
rz(-2.6601484) q[0];
rz(2.1638347) q[1];
sx q[1];
rz(-2.4974186) q[1];
sx q[1];
rz(2.0668623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9572789) q[0];
sx q[0];
rz(-1.1536351) q[0];
sx q[0];
rz(0.29846141) q[0];
rz(-pi) q[1];
rz(0.75047173) q[2];
sx q[2];
rz(-0.40130645) q[2];
sx q[2];
rz(-1.4913781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8930298) q[1];
sx q[1];
rz(-2.6932635) q[1];
sx q[1];
rz(2.1488229) q[1];
rz(0.67258622) q[3];
sx q[3];
rz(-2.6574316) q[3];
sx q[3];
rz(-1.8046093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1017199) q[2];
sx q[2];
rz(-1.1168672) q[2];
sx q[2];
rz(-2.5679585) q[2];
rz(-1.6576069) q[3];
sx q[3];
rz(-1.0480169) q[3];
sx q[3];
rz(2.535533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1178591) q[0];
sx q[0];
rz(-0.25256279) q[0];
sx q[0];
rz(-0.31164393) q[0];
rz(-0.32870865) q[1];
sx q[1];
rz(-1.3361822) q[1];
sx q[1];
rz(0.74105826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6481095) q[0];
sx q[0];
rz(-2.3435763) q[0];
sx q[0];
rz(0.9106967) q[0];
rz(-1.8297878) q[2];
sx q[2];
rz(-2.9051068) q[2];
sx q[2];
rz(1.2843997) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5697244) q[1];
sx q[1];
rz(-2.2510438) q[1];
sx q[1];
rz(0.081357439) q[1];
rz(-pi) q[2];
rz(3.1407194) q[3];
sx q[3];
rz(-0.93993087) q[3];
sx q[3];
rz(-1.6959977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7753503) q[2];
sx q[2];
rz(-0.58997184) q[2];
sx q[2];
rz(0.74014202) q[2];
rz(-2.7548693) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(-2.6630785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1665523) q[0];
sx q[0];
rz(-0.98169011) q[0];
sx q[0];
rz(-2.8983086) q[0];
rz(-2.2443306) q[1];
sx q[1];
rz(-1.5586531) q[1];
sx q[1];
rz(0.74403393) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720111) q[0];
sx q[0];
rz(-1.9599304) q[0];
sx q[0];
rz(1.6604108) q[0];
rz(-2.0646413) q[2];
sx q[2];
rz(-2.4429446) q[2];
sx q[2];
rz(0.71997627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.084483844) q[1];
sx q[1];
rz(-1.9388102) q[1];
sx q[1];
rz(2.1294566) q[1];
x q[2];
rz(1.2598557) q[3];
sx q[3];
rz(-2.5937383) q[3];
sx q[3];
rz(-2.1825298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1870785) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(-3.1336866) q[2];
rz(-1.6451969) q[3];
sx q[3];
rz(-3.0683066) q[3];
sx q[3];
rz(2.6325398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1169443) q[0];
sx q[0];
rz(-0.032289676) q[0];
sx q[0];
rz(0.13667983) q[0];
rz(-1.7928803) q[1];
sx q[1];
rz(-1.8486479) q[1];
sx q[1];
rz(2.6332556) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240468) q[0];
sx q[0];
rz(-1.4410748) q[0];
sx q[0];
rz(-1.172439) q[0];
rz(-pi) q[1];
rz(-0.29859297) q[2];
sx q[2];
rz(-0.91690874) q[2];
sx q[2];
rz(2.8841126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66022422) q[1];
sx q[1];
rz(-1.9268039) q[1];
sx q[1];
rz(-0.86467177) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.062183) q[3];
sx q[3];
rz(-1.4520922) q[3];
sx q[3];
rz(1.0371003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4623744) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(-0.057961658) q[2];
rz(0.95311779) q[3];
sx q[3];
rz(-2.0942196) q[3];
sx q[3];
rz(-0.63547772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5304831) q[0];
sx q[0];
rz(-1.7246752) q[0];
sx q[0];
rz(0.86607754) q[0];
rz(1.7762314) q[1];
sx q[1];
rz(-1.583464) q[1];
sx q[1];
rz(0.80397111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3063076) q[0];
sx q[0];
rz(-2.270335) q[0];
sx q[0];
rz(-2.7650096) q[0];
rz(-0.44064327) q[2];
sx q[2];
rz(-2.302495) q[2];
sx q[2];
rz(-0.64873141) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.436299) q[1];
sx q[1];
rz(-1.9914522) q[1];
sx q[1];
rz(2.3068025) q[1];
x q[2];
rz(-1.1671806) q[3];
sx q[3];
rz(-1.5909373) q[3];
sx q[3];
rz(0.27475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0552401) q[2];
sx q[2];
rz(-0.1487727) q[2];
sx q[2];
rz(1.0521592) q[2];
rz(0.85938984) q[3];
sx q[3];
rz(-0.79901564) q[3];
sx q[3];
rz(-0.94256443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4683485) q[0];
sx q[0];
rz(-2.4788661) q[0];
sx q[0];
rz(-2.7224139) q[0];
rz(3.1255417) q[1];
sx q[1];
rz(-1.586986) q[1];
sx q[1];
rz(0.12414653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0932179) q[0];
sx q[0];
rz(-1.7233053) q[0];
sx q[0];
rz(1.0145864) q[0];
rz(-pi) q[1];
rz(-0.89416482) q[2];
sx q[2];
rz(-2.8196206) q[2];
sx q[2];
rz(-1.3077259) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.032822306) q[1];
sx q[1];
rz(-1.3941996) q[1];
sx q[1];
rz(-2.887421) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8087555) q[3];
sx q[3];
rz(-0.66413022) q[3];
sx q[3];
rz(0.40011621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.19196709) q[2];
sx q[2];
rz(-0.73835915) q[2];
sx q[2];
rz(-0.16853608) q[2];
rz(-1.4847697) q[3];
sx q[3];
rz(-2.6717581) q[3];
sx q[3];
rz(0.43009871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94854245) q[0];
sx q[0];
rz(-2.6621303) q[0];
sx q[0];
rz(2.9051176) q[0];
rz(0.82462689) q[1];
sx q[1];
rz(-1.6025447) q[1];
sx q[1];
rz(-1.2784169) q[1];
rz(-2.5599418) q[2];
sx q[2];
rz(-0.86514513) q[2];
sx q[2];
rz(-1.9164597) q[2];
rz(-1.6705728) q[3];
sx q[3];
rz(-2.592516) q[3];
sx q[3];
rz(1.9883131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
