OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.934259) q[0];
sx q[0];
rz(-0.59036314) q[0];
sx q[0];
rz(-2.7705749) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24297548) q[0];
sx q[0];
rz(-0.78865047) q[0];
sx q[0];
rz(0.9289766) q[0];
x q[1];
rz(-2.0618093) q[2];
sx q[2];
rz(-2.2248189) q[2];
sx q[2];
rz(-0.71066463) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10378097) q[1];
sx q[1];
rz(-1.722464) q[1];
sx q[1];
rz(-1.9019466) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0891221) q[3];
sx q[3];
rz(-0.82004181) q[3];
sx q[3];
rz(-0.47580556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(0.98891813) q[2];
rz(-0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724021) q[0];
sx q[0];
rz(-0.11226421) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.8583863) q[1];
sx q[1];
rz(-2.4172799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(2.2749167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.486859) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(2.4059911) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69246768) q[1];
sx q[1];
rz(-1.7763224) q[1];
sx q[1];
rz(1.9301901) q[1];
x q[2];
rz(-1.768126) q[3];
sx q[3];
rz(-1.337968) q[3];
sx q[3];
rz(-0.80430921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(-2.779707) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1419462) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.746159) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-2.7170083) q[1];
sx q[1];
rz(-1.2190855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31419471) q[0];
sx q[0];
rz(-0.44566804) q[0];
sx q[0];
rz(2.218194) q[0];
x q[1];
rz(0.8692603) q[2];
sx q[2];
rz(-1.7995036) q[2];
sx q[2];
rz(-1.9321835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8428426) q[1];
sx q[1];
rz(-1.0423653) q[1];
sx q[1];
rz(-0.28038402) q[1];
rz(-pi) q[2];
rz(1.3027906) q[3];
sx q[3];
rz(-1.5558814) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(2.5727663) q[3];
sx q[3];
rz(-0.84972644) q[3];
sx q[3];
rz(0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(2.6233327) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-2.3148361) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137705) q[0];
sx q[0];
rz(-2.4670521) q[0];
sx q[0];
rz(0.69112372) q[0];
x q[1];
rz(-0.13917285) q[2];
sx q[2];
rz(-0.77251245) q[2];
sx q[2];
rz(0.13538361) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0485059) q[1];
sx q[1];
rz(-0.48492453) q[1];
sx q[1];
rz(0.64694689) q[1];
rz(-pi) q[2];
rz(1.2590253) q[3];
sx q[3];
rz(-0.87083737) q[3];
sx q[3];
rz(1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9892019) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(1.3789122) q[2];
rz(-3.0692696) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(-1.4962083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82692659) q[0];
sx q[0];
rz(-0.62012726) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(0.28516969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681686) q[0];
sx q[0];
rz(-1.0963206) q[0];
sx q[0];
rz(3.0939177) q[0];
x q[1];
rz(0.6638078) q[2];
sx q[2];
rz(-1.2592578) q[2];
sx q[2];
rz(0.89154348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3030745) q[1];
sx q[1];
rz(-0.31305602) q[1];
sx q[1];
rz(1.5694373) q[1];
rz(-pi) q[2];
rz(0.22484803) q[3];
sx q[3];
rz(-2.7307011) q[3];
sx q[3];
rz(-1.9083244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8482762) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(0.53945333) q[2];
rz(2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25813112) q[0];
sx q[0];
rz(-0.75043172) q[0];
sx q[0];
rz(-2.9845797) q[0];
rz(-0.69333386) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(-1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088928662) q[0];
sx q[0];
rz(-3.0658709) q[0];
sx q[0];
rz(-1.1195539) q[0];
rz(-pi) q[1];
x q[1];
rz(0.579367) q[2];
sx q[2];
rz(-2.2420792) q[2];
sx q[2];
rz(2.6210149) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4945592) q[1];
sx q[1];
rz(-1.4899947) q[1];
sx q[1];
rz(-1.8606436) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7886224) q[3];
sx q[3];
rz(-1.6680696) q[3];
sx q[3];
rz(-2.4091165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8994609) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-2.2677299) q[0];
rz(0.44772398) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0872333) q[0];
sx q[0];
rz(-1.5938252) q[0];
sx q[0];
rz(1.5646294) q[0];
rz(-pi) q[1];
rz(-2.7295693) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(-2.3840981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9782941) q[1];
sx q[1];
rz(-0.98493176) q[1];
sx q[1];
rz(-2.3818124) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53415926) q[3];
sx q[3];
rz(-1.8338406) q[3];
sx q[3];
rz(2.0260889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0447023) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-1.0905617) q[3];
sx q[3];
rz(-2.2138514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8996745) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(0.29712594) q[0];
rz(-1.3946474) q[1];
sx q[1];
rz(-1.1509832) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32983366) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(1.1128725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4722733) q[2];
sx q[2];
rz(-1.4204645) q[2];
sx q[2];
rz(2.3538102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7568126) q[1];
sx q[1];
rz(-1.9128748) q[1];
sx q[1];
rz(2.6095819) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5313247) q[3];
sx q[3];
rz(-2.2329997) q[3];
sx q[3];
rz(0.83412795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-0.45483744) q[2];
rz(2.440195) q[3];
sx q[3];
rz(-1.0304136) q[3];
sx q[3];
rz(-1.1340244) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-13*pi/16) q[0];
sx q[0];
rz(2.3440857) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-0.8126173) q[1];
sx q[1];
rz(3.033175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-0.9066559) q[0];
x q[1];
rz(1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(-2.2968963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0260967) q[1];
sx q[1];
rz(-0.27448389) q[1];
sx q[1];
rz(0.91067578) q[1];
rz(-pi) q[2];
rz(-0.18687825) q[3];
sx q[3];
rz(-0.87516057) q[3];
sx q[3];
rz(-2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(-2.8016395) q[2];
rz(0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(-3.0864339) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.508651) q[0];
sx q[0];
rz(-1.5061437) q[0];
sx q[0];
rz(-2.7642194) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55969413) q[2];
sx q[2];
rz(-1.392138) q[2];
sx q[2];
rz(0.44911227) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.89796472) q[1];
sx q[1];
rz(-2.1568255) q[1];
sx q[1];
rz(2.2378053) q[1];
rz(-pi) q[2];
rz(0.03667128) q[3];
sx q[3];
rz(-1.1150556) q[3];
sx q[3];
rz(-2.7367221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(-2.6514163) q[2];
rz(-3.0040719) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-0.93808758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162083) q[0];
sx q[0];
rz(-1.890247) q[0];
sx q[0];
rz(2.4631137) q[0];
rz(-0.2086808) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(0.31221496) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(-0.24003868) q[3];
sx q[3];
rz(-1.6143027) q[3];
sx q[3];
rz(-1.3054813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
