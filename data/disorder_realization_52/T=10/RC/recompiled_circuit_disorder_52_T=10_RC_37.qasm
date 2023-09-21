OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(4.0951401) q[1];
sx q[1];
rz(4.6440754) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8774672) q[0];
sx q[0];
rz(-1.929951) q[0];
sx q[0];
rz(-0.19192625) q[0];
rz(-pi) q[1];
rz(0.0083562334) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(-1.9985808) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68617649) q[1];
sx q[1];
rz(-1.505291) q[1];
sx q[1];
rz(0.24645611) q[1];
rz(-pi) q[2];
rz(1.9804269) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(2.9890271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(0.65594977) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-2.7080652) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(3.1343592) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4097071) q[0];
sx q[0];
rz(-1.1792372) q[0];
sx q[0];
rz(0.10086985) q[0];
rz(-pi) q[1];
rz(0.8231926) q[2];
sx q[2];
rz(-0.41831145) q[2];
sx q[2];
rz(-0.38564607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.05939535) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(-2.5897964) q[1];
x q[2];
rz(-0.98874436) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(0.83272979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(-0.69765222) q[2];
rz(0.12156045) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(-0.26161584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25900349) q[0];
sx q[0];
rz(-0.020375229) q[0];
sx q[0];
rz(-1.0783608) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82654731) q[2];
sx q[2];
rz(-1.6625704) q[2];
sx q[2];
rz(1.9020136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2178104) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(-1.4008821) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0858585) q[3];
sx q[3];
rz(-2.627943) q[3];
sx q[3];
rz(2.0897712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(2.938081) q[2];
rz(-0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(1.571636) q[0];
rz(-1.0034026) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(1.8932231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6682537) q[0];
sx q[0];
rz(-0.11675294) q[0];
sx q[0];
rz(0.97572414) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4291971) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(-1.4272387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4444256) q[1];
sx q[1];
rz(-2.4788692) q[1];
sx q[1];
rz(-2.9602833) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8321886) q[3];
sx q[3];
rz(-0.86463118) q[3];
sx q[3];
rz(2.1690926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(0.90243375) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9426614) q[0];
sx q[0];
rz(-2.5006602) q[0];
sx q[0];
rz(-3.065227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47307737) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(-2.9830473) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24689281) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(2.488399) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2540713) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(1.4590013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(1.6852089) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-0.20176372) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1267032) q[0];
sx q[0];
rz(-0.12173437) q[0];
sx q[0];
rz(-0.99862167) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47735881) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(-1.8237643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7355431) q[1];
sx q[1];
rz(-2.4821401) q[1];
sx q[1];
rz(-2.4156648) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73268907) q[3];
sx q[3];
rz(-1.9273888) q[3];
sx q[3];
rz(0.69571146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(2.8743437) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(-0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38521117) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(0.60473196) q[0];
x q[1];
rz(2.1812181) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-2.9396217) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6616933) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(0.9049306) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0833756) q[3];
sx q[3];
rz(-1.3430542) q[3];
sx q[3];
rz(-0.77373576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1192347) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(-2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(-2.5121636) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-2.8576635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0547202) q[0];
sx q[0];
rz(-2.3754639) q[0];
sx q[0];
rz(-0.92673577) q[0];
x q[1];
rz(0.21247272) q[2];
sx q[2];
rz(-0.61056911) q[2];
sx q[2];
rz(3.0688822) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64360196) q[1];
sx q[1];
rz(-2.0703348) q[1];
sx q[1];
rz(-1.040578) q[1];
rz(1.6007401) q[3];
sx q[3];
rz(-2.2303914) q[3];
sx q[3];
rz(2.3823882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.2109057) q[2];
rz(-0.25990137) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(3.0652086) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(0.22924766) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.680826) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5979249) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(0.062505917) q[0];
x q[1];
rz(0.37947189) q[2];
sx q[2];
rz(-1.737829) q[2];
sx q[2];
rz(-0.52969474) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8742121) q[1];
sx q[1];
rz(-2.0955718) q[1];
sx q[1];
rz(1.5966148) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9577515) q[3];
sx q[3];
rz(-0.82434067) q[3];
sx q[3];
rz(-1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-2.4460068) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(-0.38063231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99285179) q[0];
sx q[0];
rz(-1.7721575) q[0];
sx q[0];
rz(-0.071417884) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98999087) q[2];
sx q[2];
rz(-0.91772807) q[2];
sx q[2];
rz(-2.4534015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8457348) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(-0.22175281) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1553414) q[3];
sx q[3];
rz(-2.1086153) q[3];
sx q[3];
rz(-2.1842063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0011065817) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(2.005119) q[2];
rz(-0.040955695) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(2.1144755) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(1.8006051) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(1.1313664) q[3];
sx q[3];
rz(-1.7527179) q[3];
sx q[3];
rz(-2.6408165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];