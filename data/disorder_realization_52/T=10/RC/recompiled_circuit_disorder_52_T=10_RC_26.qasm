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
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9031154) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(-1.9360916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56514481) q[2];
sx q[2];
rz(-1.5752715) q[2];
sx q[2];
rz(-0.43484136) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.002418) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(-1.1611657) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(2.9890271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(-2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(-2.9128089) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418966) q[0];
sx q[0];
rz(-1.6640088) q[0];
sx q[0];
rz(1.1774363) q[0];
x q[1];
rz(-2.8480808) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(-2.7345865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65456395) q[1];
sx q[1];
rz(-1.2083168) q[1];
sx q[1];
rz(-1.3351721) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90356566) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(0.69765222) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.7720222) q[0];
rz(-1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(2.8799768) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25900349) q[0];
sx q[0];
rz(-3.1212174) q[0];
sx q[0];
rz(2.0632319) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3150453) q[2];
sx q[2];
rz(-1.6625704) q[2];
sx q[2];
rz(-1.9020136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2178104) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(1.4008821) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6022127) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-2.2198548) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(2.8620499) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(1.5699566) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.8932231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49444775) q[0];
sx q[0];
rz(-1.6361423) q[0];
sx q[0];
rz(1.6676184) q[0];
x q[1];
rz(-0.7123956) q[2];
sx q[2];
rz(-2.8678896) q[2];
sx q[2];
rz(1.4272387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.6819685) q[1];
sx q[1];
rz(-2.486869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72327153) q[3];
sx q[3];
rz(-1.7687106) q[3];
sx q[3];
rz(2.7151782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(2.3045585) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(-2.3663882) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(0.90243375) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2941023) q[0];
sx q[0];
rz(-2.2095564) q[0];
sx q[0];
rz(-1.5139447) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47307737) q[2];
sx q[2];
rz(-2.6396751) q[2];
sx q[2];
rz(0.15854533) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.21197) q[1];
sx q[1];
rz(-2.3096482) q[1];
sx q[1];
rz(0.79939876) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5621454) q[3];
sx q[3];
rz(-0.31673613) q[3];
sx q[3];
rz(2.1474311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-2.6384171) q[0];
rz(1.6852089) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(2.9398289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87529463) q[0];
sx q[0];
rz(-1.6365956) q[0];
sx q[0];
rz(1.6732897) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6642338) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(-1.3178283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6946053) q[1];
sx q[1];
rz(-1.9896549) q[1];
sx q[1];
rz(2.6161731) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1062578) q[3];
sx q[3];
rz(-2.2482276) q[3];
sx q[3];
rz(2.5708452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(-2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-0.057549495) q[0];
rz(1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7563815) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(-0.60473196) q[0];
rz(-2.9808688) q[2];
sx q[2];
rz(-1.7956927) q[2];
sx q[2];
rz(2.7149534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-0.53240314) q[1];
sx q[1];
rz(0.9049306) q[1];
x q[2];
rz(-2.8816678) q[3];
sx q[3];
rz(-1.0726895) q[3];
sx q[3];
rz(0.92344027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(-2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(0.27012816) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-2.8576635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1536381) q[0];
sx q[0];
rz(-2.0001912) q[0];
sx q[0];
rz(0.9149787) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9291199) q[2];
sx q[2];
rz(-2.5310235) q[2];
sx q[2];
rz(3.0688822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9406416) q[1];
sx q[1];
rz(-2.0307699) q[1];
sx q[1];
rz(2.5775787) q[1];
rz(-pi) q[2];
rz(-1.5408526) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(0.75920446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.930687) q[2];
rz(0.25990137) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(-2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(2.912345) q[0];
rz(-2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(-1.4607666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5979249) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(-0.062505917) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3912348) q[2];
sx q[2];
rz(-1.1968687) q[2];
sx q[2];
rz(-1.1073081) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9257099) q[1];
sx q[1];
rz(-0.52535086) q[1];
sx q[1];
rz(3.0970296) q[1];
rz(0.38735729) q[3];
sx q[3];
rz(-2.318317) q[3];
sx q[3];
rz(0.92126095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-2.4460068) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(-2.8046872) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-0.38063231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363627) q[0];
sx q[0];
rz(-0.21348937) q[0];
sx q[0];
rz(1.9070894) q[0];
rz(-pi) q[1];
rz(0.98999087) q[2];
sx q[2];
rz(-2.2238646) q[2];
sx q[2];
rz(-2.4534015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29585782) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(0.22175281) q[1];
rz(-pi) q[2];
rz(-1.9862513) q[3];
sx q[3];
rz(-2.1086153) q[3];
sx q[3];
rz(-2.1842063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(2.005119) q[2];
rz(-0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.3142746) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(-0.2480416) q[2];
sx q[2];
rz(-1.7938062) q[2];
sx q[2];
rz(1.2563406) q[2];
rz(1.9789226) q[3];
sx q[3];
rz(-0.47331953) q[3];
sx q[3];
rz(2.4389653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
