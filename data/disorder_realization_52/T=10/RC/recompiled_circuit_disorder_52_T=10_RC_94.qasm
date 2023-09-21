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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8774672) q[0];
sx q[0];
rz(-1.929951) q[0];
sx q[0];
rz(2.9496664) q[0];
rz(-pi) q[1];
rz(-0.56514481) q[2];
sx q[2];
rz(-1.5663212) q[2];
sx q[2];
rz(-2.7067513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1391746) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7847071) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(-0.68912904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(-0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475964) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418966) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.9641563) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29351182) q[2];
sx q[2];
rz(-1.268317) q[2];
sx q[2];
rz(2.7345865) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0821973) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(0.55179623) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1528483) q[3];
sx q[3];
rz(-0.75609222) q[3];
sx q[3];
rz(-0.83272979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(0.69765222) q[2];
rz(-0.12156045) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(-2.8799768) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944377) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(-1.5887512) q[0];
rz(1.4357655) q[2];
sx q[2];
rz(-0.74880744) q[2];
sx q[2];
rz(2.7111862) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9237823) q[1];
sx q[1];
rz(-2.3465956) q[1];
sx q[1];
rz(1.4008821) q[1];
rz(-pi) q[2];
rz(-0.055734169) q[3];
sx q[3];
rz(-2.627943) q[3];
sx q[3];
rz(1.0518215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-2.2198548) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.8932231) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715863) q[0];
sx q[0];
rz(-1.6674111) q[0];
sx q[0];
rz(-0.065652547) q[0];
rz(2.4291971) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(-1.4272387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9256546) q[1];
sx q[1];
rz(-2.2207894) q[1];
sx q[1];
rz(-1.4309806) q[1];
x q[2];
rz(-0.29420935) q[3];
sx q[3];
rz(-2.3964786) q[3];
sx q[3];
rz(-1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0288329) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0060624881) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(-0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9426614) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(-3.065227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4544223) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(0.99046747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8946998) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(2.488399) q[1];
rz(3.1387572) q[3];
sx q[3];
rz(-1.8875202) q[3];
sx q[3];
rz(1.0032652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9465785) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(-1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.4590013) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3301795) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-1.0738942) q[1];
sx q[1];
rz(2.9398289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70226442) q[0];
sx q[0];
rz(-1.6730671) q[0];
sx q[0];
rz(3.0754473) q[0];
rz(-pi) q[1];
rz(2.0189507) q[2];
sx q[2];
rz(-2.0070224) q[2];
sx q[2];
rz(-3.089038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.7355431) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(-0.72592782) q[1];
x q[2];
rz(2.4089036) q[3];
sx q[3];
rz(-1.2142039) q[3];
sx q[3];
rz(-2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(-2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(3.0840432) q[0];
rz(-1.6607704) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38521117) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(0.60473196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16072388) q[2];
sx q[2];
rz(-1.3458999) q[2];
sx q[2];
rz(-2.7149534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2191636) q[1];
sx q[1];
rz(-1.1601829) q[1];
sx q[1];
rz(-2.7925522) q[1];
x q[2];
rz(-1.1293291) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(-0.41551513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-2.4286178) q[1];
sx q[1];
rz(2.8576635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89284183) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(-0.52389223) q[0];
x q[1];
rz(2.9291199) q[2];
sx q[2];
rz(-0.61056911) q[2];
sx q[2];
rz(-3.0688822) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4979907) q[1];
sx q[1];
rz(-1.0712578) q[1];
sx q[1];
rz(1.040578) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-2.4814197) q[3];
sx q[3];
rz(2.431228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55390629) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(1.930687) q[2];
rz(-0.25990137) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54366771) q[0];
sx q[0];
rz(-2.857508) q[0];
sx q[0];
rz(0.062505917) q[0];
x q[1];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(-1.1073081) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31635346) q[1];
sx q[1];
rz(-1.5931399) q[1];
sx q[1];
rz(0.52492001) q[1];
rz(1.9577515) q[3];
sx q[3];
rz(-0.82434067) q[3];
sx q[3];
rz(1.4617621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(0.6955859) q[2];
rz(-0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5678976) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(-2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(2.7609603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99285179) q[0];
sx q[0];
rz(-1.3694351) q[0];
sx q[0];
rz(0.071417884) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5194174) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(-0.13571339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29585782) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(2.9198398) q[1];
x q[2];
rz(0.57772824) q[3];
sx q[3];
rz(-1.9247705) q[3];
sx q[3];
rz(2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(2.005119) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363591) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
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
rz(-1.1313664) q[3];
sx q[3];
rz(-1.3888748) q[3];
sx q[3];
rz(0.5007762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
