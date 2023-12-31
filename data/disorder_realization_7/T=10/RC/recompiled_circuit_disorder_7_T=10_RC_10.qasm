OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.964736) q[0];
sx q[0];
rz(-2.2778947) q[0];
sx q[0];
rz(2.9404844) q[0];
rz(-1.3970628) q[1];
sx q[1];
rz(-1.8048598) q[1];
sx q[1];
rz(2.5147658) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39591046) q[0];
sx q[0];
rz(-0.22326176) q[0];
sx q[0];
rz(-2.1366871) q[0];
rz(-pi) q[1];
rz(-2.5457382) q[2];
sx q[2];
rz(-0.41759727) q[2];
sx q[2];
rz(-0.28996224) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1498897) q[1];
sx q[1];
rz(-0.68421364) q[1];
sx q[1];
rz(1.0690772) q[1];
rz(-pi) q[2];
rz(2.1342282) q[3];
sx q[3];
rz(-2.8197188) q[3];
sx q[3];
rz(0.25062996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8864002) q[2];
sx q[2];
rz(-2.2108086) q[2];
sx q[2];
rz(1.8480776) q[2];
rz(2.0375997) q[3];
sx q[3];
rz(-0.84775001) q[3];
sx q[3];
rz(2.3867992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8936546) q[0];
sx q[0];
rz(-1.0590483) q[0];
sx q[0];
rz(-0.98518103) q[0];
rz(2.7424116) q[1];
sx q[1];
rz(-1.8789411) q[1];
sx q[1];
rz(3.0342297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3283008) q[0];
sx q[0];
rz(-1.0257226) q[0];
sx q[0];
rz(2.6585447) q[0];
rz(-2.5537002) q[2];
sx q[2];
rz(-1.8111472) q[2];
sx q[2];
rz(2.8628778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0641891) q[1];
sx q[1];
rz(-1.2014376) q[1];
sx q[1];
rz(0.87537745) q[1];
rz(-pi) q[2];
rz(-0.3932088) q[3];
sx q[3];
rz(-2.7244096) q[3];
sx q[3];
rz(1.259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62464109) q[2];
sx q[2];
rz(-1.7525643) q[2];
sx q[2];
rz(-2.4798933) q[2];
rz(3.0858357) q[3];
sx q[3];
rz(-0.35900933) q[3];
sx q[3];
rz(-2.8957446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702883) q[0];
sx q[0];
rz(-2.5876973) q[0];
sx q[0];
rz(0.04034986) q[0];
rz(-2.5616052) q[1];
sx q[1];
rz(-2.2433387) q[1];
sx q[1];
rz(-2.0592164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4989935) q[0];
sx q[0];
rz(-1.9369619) q[0];
sx q[0];
rz(-0.55143349) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9897887) q[2];
sx q[2];
rz(-1.2200583) q[2];
sx q[2];
rz(1.1317859) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97835097) q[1];
sx q[1];
rz(-2.8132073) q[1];
sx q[1];
rz(-0.82113012) q[1];
rz(-pi) q[2];
rz(0.97614395) q[3];
sx q[3];
rz(-0.38446063) q[3];
sx q[3];
rz(0.36851766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0064156) q[2];
sx q[2];
rz(-1.2536896) q[2];
sx q[2];
rz(0.55389261) q[2];
rz(1.6484377) q[3];
sx q[3];
rz(-1.0529543) q[3];
sx q[3];
rz(-2.5890787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64788139) q[0];
sx q[0];
rz(-0.58409062) q[0];
sx q[0];
rz(-0.035382263) q[0];
rz(2.1454051) q[1];
sx q[1];
rz(-1.532998) q[1];
sx q[1];
rz(2.6534973) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6513718) q[0];
sx q[0];
rz(-2.0351366) q[0];
sx q[0];
rz(3.0547633) q[0];
x q[1];
rz(-0.57882611) q[2];
sx q[2];
rz(-1.0597214) q[2];
sx q[2];
rz(1.5103112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6285703) q[1];
sx q[1];
rz(-1.4161308) q[1];
sx q[1];
rz(-0.13048529) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.048269) q[3];
sx q[3];
rz(-2.100088) q[3];
sx q[3];
rz(0.57100163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7867243) q[2];
sx q[2];
rz(-1.0151261) q[2];
sx q[2];
rz(0.50746894) q[2];
rz(-1.9594225) q[3];
sx q[3];
rz(-1.2927262) q[3];
sx q[3];
rz(-1.2549843) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59947157) q[0];
sx q[0];
rz(-2.4276955) q[0];
sx q[0];
rz(0.3702634) q[0];
rz(-0.26369357) q[1];
sx q[1];
rz(-1.5779243) q[1];
sx q[1];
rz(-2.945074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4465966) q[0];
sx q[0];
rz(-1.8509522) q[0];
sx q[0];
rz(-0.3145991) q[0];
rz(-3.1095805) q[2];
sx q[2];
rz(-1.1099166) q[2];
sx q[2];
rz(-0.33547685) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0177808) q[1];
sx q[1];
rz(-1.8291626) q[1];
sx q[1];
rz(1.7523132) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27627857) q[3];
sx q[3];
rz(-1.7344788) q[3];
sx q[3];
rz(0.39556634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.013441) q[2];
sx q[2];
rz(-2.1898966) q[2];
sx q[2];
rz(-2.2244942) q[2];
rz(-0.83459485) q[3];
sx q[3];
rz(-1.3093427) q[3];
sx q[3];
rz(0.33368567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.82814) q[0];
sx q[0];
rz(-1.4414635) q[0];
sx q[0];
rz(0.20203461) q[0];
rz(1.0728041) q[1];
sx q[1];
rz(-1.5067116) q[1];
sx q[1];
rz(-1.7477759) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4112339) q[0];
sx q[0];
rz(-0.3404091) q[0];
sx q[0];
rz(-2.911724) q[0];
rz(-pi) q[1];
rz(-1.2932111) q[2];
sx q[2];
rz(-0.38799122) q[2];
sx q[2];
rz(-0.1462305) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.60336941) q[1];
sx q[1];
rz(-1.9586316) q[1];
sx q[1];
rz(-2.8549854) q[1];
rz(0.73563852) q[3];
sx q[3];
rz(-1.9197575) q[3];
sx q[3];
rz(-0.58128202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1918734) q[2];
sx q[2];
rz(-1.4800025) q[2];
sx q[2];
rz(-2.2992772) q[2];
rz(2.0541644) q[3];
sx q[3];
rz(-2.4098586) q[3];
sx q[3];
rz(1.6585763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52952805) q[0];
sx q[0];
rz(-1.8719215) q[0];
sx q[0];
rz(-2.1202309) q[0];
rz(1.1313324) q[1];
sx q[1];
rz(-0.94968692) q[1];
sx q[1];
rz(-2.9313415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049858) q[0];
sx q[0];
rz(-0.14963089) q[0];
sx q[0];
rz(2.0307226) q[0];
rz(0.37961752) q[2];
sx q[2];
rz(-2.9656177) q[2];
sx q[2];
rz(1.9642252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6107632) q[1];
sx q[1];
rz(-1.7621343) q[1];
sx q[1];
rz(-2.8192239) q[1];
rz(-pi) q[2];
rz(-2.290091) q[3];
sx q[3];
rz(-1.1151033) q[3];
sx q[3];
rz(1.15629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3147099) q[2];
sx q[2];
rz(-1.751739) q[2];
sx q[2];
rz(-1.7604527) q[2];
rz(0.76210493) q[3];
sx q[3];
rz(-0.47587454) q[3];
sx q[3];
rz(2.7281318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49522266) q[0];
sx q[0];
rz(-2.3587527) q[0];
sx q[0];
rz(0.58724171) q[0];
rz(0.44627055) q[1];
sx q[1];
rz(-1.279436) q[1];
sx q[1];
rz(2.1648724) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2864838) q[0];
sx q[0];
rz(-2.6628462) q[0];
sx q[0];
rz(-0.32740645) q[0];
rz(-pi) q[1];
rz(-3.061053) q[2];
sx q[2];
rz(-1.7320247) q[2];
sx q[2];
rz(1.5494407) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50423065) q[1];
sx q[1];
rz(-1.1334051) q[1];
sx q[1];
rz(0.38423844) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5891853) q[3];
sx q[3];
rz(-2.2033785) q[3];
sx q[3];
rz(1.4668902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7405159) q[2];
sx q[2];
rz(-2.433321) q[2];
sx q[2];
rz(-2.5891499) q[2];
rz(-0.46594122) q[3];
sx q[3];
rz(-1.3876785) q[3];
sx q[3];
rz(2.1140816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4733646) q[0];
sx q[0];
rz(-2.3478274) q[0];
sx q[0];
rz(2.2209432) q[0];
rz(-0.051041516) q[1];
sx q[1];
rz(-1.5565245) q[1];
sx q[1];
rz(-2.2424973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63910045) q[0];
sx q[0];
rz(-0.99196767) q[0];
sx q[0];
rz(2.6274908) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0752418) q[2];
sx q[2];
rz(-0.62014183) q[2];
sx q[2];
rz(1.3401741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2832665) q[1];
sx q[1];
rz(-2.7910633) q[1];
sx q[1];
rz(1.2298898) q[1];
rz(-pi) q[2];
rz(-0.00021342834) q[3];
sx q[3];
rz(-1.9409688) q[3];
sx q[3];
rz(-2.7406524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1961394) q[2];
sx q[2];
rz(-3.0710199) q[2];
sx q[2];
rz(3.0370039) q[2];
rz(0.83834046) q[3];
sx q[3];
rz(-0.74931562) q[3];
sx q[3];
rz(-2.2588363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8024682) q[0];
sx q[0];
rz(-2.323928) q[0];
sx q[0];
rz(1.1178281) q[0];
rz(2.4291908) q[1];
sx q[1];
rz(-0.63122216) q[1];
sx q[1];
rz(-0.39696473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8660276) q[0];
sx q[0];
rz(-2.4336928) q[0];
sx q[0];
rz(2.3562354) q[0];
rz(-pi) q[1];
rz(2.2592779) q[2];
sx q[2];
rz(-0.46122641) q[2];
sx q[2];
rz(-0.099909401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1399999) q[1];
sx q[1];
rz(-1.6591751) q[1];
sx q[1];
rz(-1.9077076) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40542094) q[3];
sx q[3];
rz(-1.1076895) q[3];
sx q[3];
rz(-2.6904358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1931856) q[2];
sx q[2];
rz(-2.2414175) q[2];
sx q[2];
rz(-0.096654264) q[2];
rz(-1.4139253) q[3];
sx q[3];
rz(-2.0035412) q[3];
sx q[3];
rz(-2.0146501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5554572) q[0];
sx q[0];
rz(-2.0270892) q[0];
sx q[0];
rz(-2.8591697) q[0];
rz(1.8718406) q[1];
sx q[1];
rz(-1.3602263) q[1];
sx q[1];
rz(0.7855986) q[1];
rz(1.6554228) q[2];
sx q[2];
rz(-1.6783236) q[2];
sx q[2];
rz(2.8477737) q[2];
rz(-1.312064) q[3];
sx q[3];
rz(-1.6862292) q[3];
sx q[3];
rz(-0.31173691) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
