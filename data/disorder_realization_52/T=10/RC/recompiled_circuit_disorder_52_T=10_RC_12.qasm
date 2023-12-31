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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2384773) q[0];
sx q[0];
rz(-1.7503386) q[0];
sx q[0];
rz(-1.205501) q[0];
rz(-1.5654972) q[2];
sx q[2];
rz(-2.1359348) q[2];
sx q[2];
rz(1.1331171) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.002418) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-0.26267085) q[1];
rz(-pi) q[2];
rz(-1.9804269) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(-2.9890271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475964) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(-2.7080652) q[0];
rz(0.22878376) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(-0.0072335009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418966) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.1774363) q[0];
rz(-pi) q[1];
rz(-0.29351182) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(-0.40700618) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65456395) q[1];
sx q[1];
rz(-1.2083168) q[1];
sx q[1];
rz(-1.3351721) q[1];
rz(-pi) q[2];
rz(-0.47827999) q[3];
sx q[3];
rz(-0.96049958) q[3];
sx q[3];
rz(1.573521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5269512) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6737297) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(1.7720222) q[0];
rz(-1.2415775) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(0.26161584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25900349) q[0];
sx q[0];
rz(-0.020375229) q[0];
sx q[0];
rz(-2.0632319) q[0];
x q[1];
rz(0.82654731) q[2];
sx q[2];
rz(-1.4790223) q[2];
sx q[2];
rz(1.2395791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9237823) q[1];
sx q[1];
rz(-0.79499704) q[1];
sx q[1];
rz(-1.7407106) q[1];
rz(0.055734169) q[3];
sx q[3];
rz(-2.627943) q[3];
sx q[3];
rz(2.0897712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(2.2198548) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(1.2483695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49444775) q[0];
sx q[0];
rz(-1.5054504) q[0];
sx q[0];
rz(-1.6676184) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7123956) q[2];
sx q[2];
rz(-2.8678896) q[2];
sx q[2];
rz(1.7143539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9256546) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(1.4309806) q[1];
rz(0.72327153) q[3];
sx q[3];
rz(-1.372882) q[3];
sx q[3];
rz(-0.42641446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9426614) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(-3.065227) q[0];
x q[1];
rz(2.6871704) q[2];
sx q[2];
rz(-1.3497958) q[2];
sx q[2];
rz(2.1511252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9296226) q[1];
sx q[1];
rz(-0.8319444) q[1];
sx q[1];
rz(-2.3421939) q[1];
rz(0.0028354672) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(2.1383274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.1966594) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3301795) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(0.50317558) q[0];
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
rz(-0.87529463) q[0];
sx q[0];
rz(-1.504997) q[0];
sx q[0];
rz(-1.6732897) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0189507) q[2];
sx q[2];
rz(-1.1345703) q[2];
sx q[2];
rz(-3.089038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10776785) q[1];
sx q[1];
rz(-2.0467842) q[1];
sx q[1];
rz(-1.0955217) q[1];
rz(-pi) q[2];
rz(2.6334409) q[3];
sx q[3];
rz(-0.80012946) q[3];
sx q[3];
rz(1.8964163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-3.0624793) q[3];
sx q[3];
rz(2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7563815) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-0.60473196) q[0];
rz(-2.9808688) q[2];
sx q[2];
rz(-1.7956927) q[2];
sx q[2];
rz(-0.42663923) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-0.53240314) q[1];
sx q[1];
rz(-2.2366621) q[1];
rz(2.8816678) q[3];
sx q[3];
rz(-1.0726895) q[3];
sx q[3];
rz(2.2181524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1192347) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(2.5121636) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-0.28392917) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98795456) q[0];
sx q[0];
rz(-1.1414014) q[0];
sx q[0];
rz(2.226614) q[0];
rz(1.7173041) q[2];
sx q[2];
rz(-2.1657145) q[2];
sx q[2];
rz(2.811424) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9406416) q[1];
sx q[1];
rz(-1.1108228) q[1];
sx q[1];
rz(-2.5775787) q[1];
x q[2];
rz(-3.1030032) q[3];
sx q[3];
rz(-2.4814197) q[3];
sx q[3];
rz(0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.4607666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0544573) q[0];
sx q[0];
rz(-1.5883049) q[0];
sx q[0];
rz(-2.8580335) q[0];
rz(-0.42713366) q[2];
sx q[2];
rz(-0.41296994) q[2];
sx q[2];
rz(0.6461179) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9257099) q[1];
sx q[1];
rz(-0.52535086) q[1];
sx q[1];
rz(-3.0970296) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78482307) q[3];
sx q[3];
rz(-1.8514957) q[3];
sx q[3];
rz(-0.37898889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(2.4263884) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.573695) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(-2.8046872) q[0];
rz(-0.20740549) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1487409) q[0];
sx q[0];
rz(-1.7721575) q[0];
sx q[0];
rz(-0.071417884) q[0];
x q[1];
rz(0.98999087) q[2];
sx q[2];
rz(-2.2238646) q[2];
sx q[2];
rz(-2.4534015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29585782) q[1];
sx q[1];
rz(-2.0669792) q[1];
sx q[1];
rz(2.9198398) q[1];
rz(-pi) q[2];
rz(1.1553414) q[3];
sx q[3];
rz(-2.1086153) q[3];
sx q[3];
rz(0.95738639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(1.1364737) q[2];
rz(-0.040955695) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052335) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(-1.0271172) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(2.3958191) q[2];
sx q[2];
rz(-2.8095828) q[2];
sx q[2];
rz(-2.7381894) q[2];
rz(-0.20053486) q[3];
sx q[3];
rz(-1.1391098) q[3];
sx q[3];
rz(-1.1548635) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
