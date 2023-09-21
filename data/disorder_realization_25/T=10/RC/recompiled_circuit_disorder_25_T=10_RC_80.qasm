OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1782921) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(-0.51142366) q[0];
x q[1];
rz(-2.232993) q[2];
sx q[2];
rz(-2.5052921) q[2];
sx q[2];
rz(1.5696021) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3808448) q[1];
sx q[1];
rz(-1.5686085) q[1];
sx q[1];
rz(-2.6373177) q[1];
rz(1.3120996) q[3];
sx q[3];
rz(-1.5651349) q[3];
sx q[3];
rz(0.86710801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-0.56935707) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-2.5053535) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(2.5879481) q[0];
rz(-1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419932) q[0];
sx q[0];
rz(-2.2337231) q[0];
sx q[0];
rz(-2.2671188) q[0];
x q[1];
rz(0.1594752) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-2.8603539) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0765637) q[1];
sx q[1];
rz(-0.20670465) q[1];
sx q[1];
rz(-1.6480584) q[1];
rz(0.66744653) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(0.28039704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0040940293) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(2.3705204) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035633798) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(1.8133624) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0686915) q[2];
sx q[2];
rz(-2.6325668) q[2];
sx q[2];
rz(-0.44391649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33176955) q[1];
sx q[1];
rz(-0.66546813) q[1];
sx q[1];
rz(-0.21824093) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86665385) q[3];
sx q[3];
rz(-2.1979077) q[3];
sx q[3];
rz(1.31124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(2.3160034) q[0];
rz(-1.1666974) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.694214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7839514) q[0];
sx q[0];
rz(-0.28342993) q[0];
sx q[0];
rz(-2.1901603) q[0];
rz(-pi) q[1];
rz(-1.5201735) q[2];
sx q[2];
rz(-2.9630337) q[2];
sx q[2];
rz(-0.52390097) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4923258) q[1];
sx q[1];
rz(-0.71486231) q[1];
sx q[1];
rz(-1.3267172) q[1];
rz(2.7778266) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(1.0432537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.2129983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1367462) q[0];
sx q[0];
rz(-2.6803225) q[0];
sx q[0];
rz(-0.61704163) q[0];
rz(-2.5118675) q[2];
sx q[2];
rz(-2.1246398) q[2];
sx q[2];
rz(2.6385006) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0141543) q[1];
sx q[1];
rz(-1.452983) q[1];
sx q[1];
rz(-0.25728667) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9807068) q[3];
sx q[3];
rz(-2.8299035) q[3];
sx q[3];
rz(0.89524549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.1155788) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(0.58475959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11767865) q[0];
sx q[0];
rz(-2.4524134) q[0];
sx q[0];
rz(-2.6373449) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5271687) q[2];
sx q[2];
rz(-1.5922058) q[2];
sx q[2];
rz(-1.1504088) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11440052) q[1];
sx q[1];
rz(-2.235734) q[1];
sx q[1];
rz(-0.36892885) q[1];
rz(-1.5201969) q[3];
sx q[3];
rz(-1.4046602) q[3];
sx q[3];
rz(-1.7069526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(0.60639492) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43019766) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(-0.27336143) q[0];
rz(-2.0386001) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(1.4008092) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8555774) q[1];
sx q[1];
rz(-0.77076036) q[1];
sx q[1];
rz(-1.7325749) q[1];
rz(-pi) q[2];
rz(-2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(0.94999653) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-0.22769134) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3823711) q[0];
sx q[0];
rz(-1.4996487) q[0];
sx q[0];
rz(-1.9842149) q[0];
rz(-pi) q[1];
rz(1.20485) q[2];
sx q[2];
rz(-1.1507251) q[2];
sx q[2];
rz(2.7811188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5880809) q[1];
sx q[1];
rz(-1.3462102) q[1];
sx q[1];
rz(-2.5976546) q[1];
x q[2];
rz(-1.6119192) q[3];
sx q[3];
rz(-2.1293981) q[3];
sx q[3];
rz(-0.55762824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.7611327) q[3];
sx q[3];
rz(0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39695981) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(0.80378419) q[0];
rz(2.0869758) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(-1.9931591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52590695) q[0];
sx q[0];
rz(-1.6427759) q[0];
sx q[0];
rz(0.023948897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0954275) q[2];
sx q[2];
rz(-1.8881919) q[2];
sx q[2];
rz(1.5898926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0828404) q[1];
sx q[1];
rz(-1.0995004) q[1];
sx q[1];
rz(-0.48401644) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7896342) q[3];
sx q[3];
rz(-1.5378012) q[3];
sx q[3];
rz(-1.2034504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(0.7397488) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(-2.443312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88903504) q[0];
sx q[0];
rz(-1.6425942) q[0];
sx q[0];
rz(-0.30307146) q[0];
x q[1];
rz(0.024784879) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(0.63400808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2472898) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(1.1029907) q[1];
rz(1.2020338) q[3];
sx q[3];
rz(-2.5420815) q[3];
sx q[3];
rz(-2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(0.84021604) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-2.2655723) q[3];
sx q[3];
rz(1.9742112) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(3.1124658) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(1.9188948) q[2];
sx q[2];
rz(-2.1445027) q[2];
sx q[2];
rz(2.4377844) q[2];
rz(2.2610353) q[3];
sx q[3];
rz(-1.4546483) q[3];
sx q[3];
rz(1.6402257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
