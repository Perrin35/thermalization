OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5470619) q[0];
sx q[0];
rz(-2.0201595) q[0];
sx q[0];
rz(2.9601331) q[0];
rz(-1.0815066) q[1];
sx q[1];
rz(-2.4681611) q[1];
sx q[1];
rz(-2.0884617) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15381972) q[0];
sx q[0];
rz(-1.1115371) q[0];
sx q[0];
rz(0.23621969) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40197576) q[2];
sx q[2];
rz(-2.4789414) q[2];
sx q[2];
rz(1.0647917) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6932106) q[1];
sx q[1];
rz(-2.3464077) q[1];
sx q[1];
rz(0.27547142) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0041204) q[3];
sx q[3];
rz(-2.7702799) q[3];
sx q[3];
rz(0.565688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(-1.0129499) q[3];
sx q[3];
rz(-2.2949341) q[3];
sx q[3];
rz(-0.071454123) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0435836) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(1.0455421) q[1];
sx q[1];
rz(-1.7273993) q[1];
sx q[1];
rz(-1.5140623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428974) q[0];
sx q[0];
rz(-1.1637582) q[0];
sx q[0];
rz(1.0242978) q[0];
rz(-2.3081231) q[2];
sx q[2];
rz(-1.5499299) q[2];
sx q[2];
rz(-1.8881063) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82799998) q[1];
sx q[1];
rz(-2.8803821) q[1];
sx q[1];
rz(2.2952495) q[1];
rz(-0.049948143) q[3];
sx q[3];
rz(-1.4686099) q[3];
sx q[3];
rz(-3.0199043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4864768) q[2];
sx q[2];
rz(-1.0412419) q[2];
sx q[2];
rz(-2.3454323) q[2];
rz(-0.97186175) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(0.17175737) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7650836) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(-2.8308716) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.6832738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761988) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(1.7182299) q[0];
rz(2.1737353) q[2];
sx q[2];
rz(-1.0087011) q[2];
sx q[2];
rz(-2.2130655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2787784) q[1];
sx q[1];
rz(-1.9796951) q[1];
sx q[1];
rz(1.6294953) q[1];
rz(-1.5879257) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(-0.67908755) q[2];
rz(-1.7689765) q[3];
sx q[3];
rz(-1.2981828) q[3];
sx q[3];
rz(-2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(-2.0171719) q[0];
rz(1.2202948) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(1.4845928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.5201709) q[0];
sx q[0];
rz(3.1111858) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9001311) q[2];
sx q[2];
rz(-0.54154684) q[2];
sx q[2];
rz(2.3878218) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93814072) q[1];
sx q[1];
rz(-1.005299) q[1];
sx q[1];
rz(-0.52312619) q[1];
x q[2];
rz(-1.2403537) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(0.89360039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1241887) q[2];
sx q[2];
rz(-1.8576531) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(1.7381564) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6338585) q[0];
sx q[0];
rz(-1.2861179) q[0];
sx q[0];
rz(2.741709) q[0];
rz(1.1625066) q[1];
sx q[1];
rz(-1.8116654) q[1];
sx q[1];
rz(-2.9679325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0608873) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(-1.5990431) q[0];
rz(-pi) q[1];
rz(2.747614) q[2];
sx q[2];
rz(-1.7658965) q[2];
sx q[2];
rz(2.5815798) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0581936) q[1];
sx q[1];
rz(-0.10903437) q[1];
sx q[1];
rz(-1.3074247) q[1];
rz(-2.5912656) q[3];
sx q[3];
rz(-0.384207) q[3];
sx q[3];
rz(-1.3889544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6614723) q[2];
sx q[2];
rz(-1.1602594) q[2];
sx q[2];
rz(-0.76812569) q[2];
rz(-0.85401946) q[3];
sx q[3];
rz(-1.4203527) q[3];
sx q[3];
rz(0.97222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1059234) q[0];
sx q[0];
rz(-0.24374715) q[0];
sx q[0];
rz(-1.4703898) q[0];
rz(-0.51180965) q[1];
sx q[1];
rz(-0.51135951) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2238732) q[0];
sx q[0];
rz(-1.4871162) q[0];
sx q[0];
rz(-2.9907945) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6799455) q[2];
sx q[2];
rz(-0.6558154) q[2];
sx q[2];
rz(1.1193502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9547564) q[1];
sx q[1];
rz(-1.7446767) q[1];
sx q[1];
rz(-1.9435005) q[1];
x q[2];
rz(-2.0538123) q[3];
sx q[3];
rz(-2.5452168) q[3];
sx q[3];
rz(-3.0697825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64289552) q[2];
sx q[2];
rz(-2.2968447) q[2];
sx q[2];
rz(1.6112304) q[2];
rz(1.6879843) q[3];
sx q[3];
rz(-2.0724824) q[3];
sx q[3];
rz(0.24916515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11944184) q[0];
sx q[0];
rz(-0.74437737) q[0];
sx q[0];
rz(1.7013593) q[0];
rz(-2.4123736) q[1];
sx q[1];
rz(-1.9897285) q[1];
sx q[1];
rz(-1.1332606) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5232552) q[0];
sx q[0];
rz(-1.6132857) q[0];
sx q[0];
rz(2.8190814) q[0];
rz(-pi) q[1];
rz(-2.4293443) q[2];
sx q[2];
rz(-1.5603666) q[2];
sx q[2];
rz(1.041009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3164191) q[1];
sx q[1];
rz(-0.9269886) q[1];
sx q[1];
rz(-0.34367798) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6071068) q[3];
sx q[3];
rz(-2.839698) q[3];
sx q[3];
rz(-1.2045977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7363654) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-0.48842946) q[2];
rz(-1.8296261) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(-0.64129889) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(1.2088998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9933388) q[0];
sx q[0];
rz(-0.78972602) q[0];
sx q[0];
rz(-1.2475616) q[0];
x q[1];
rz(0.18431633) q[2];
sx q[2];
rz(-1.3375998) q[2];
sx q[2];
rz(-1.2554982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88747957) q[1];
sx q[1];
rz(-1.5548318) q[1];
sx q[1];
rz(-1.70114) q[1];
rz(-pi) q[2];
rz(1.204793) q[3];
sx q[3];
rz(-0.98099785) q[3];
sx q[3];
rz(0.78192151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(1.570638) q[2];
rz(-2.2682244) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-2.7895555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97380012) q[0];
sx q[0];
rz(-1.5252824) q[0];
sx q[0];
rz(-2.8299676) q[0];
rz(-2.3198126) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.5100381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8845997) q[0];
sx q[0];
rz(-1.8879226) q[0];
sx q[0];
rz(1.3604128) q[0];
rz(-pi) q[1];
rz(-1.0636343) q[2];
sx q[2];
rz(-2.7535451) q[2];
sx q[2];
rz(-1.0840814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1875302) q[1];
sx q[1];
rz(-0.64779753) q[1];
sx q[1];
rz(2.8026583) q[1];
rz(-pi) q[2];
rz(1.7543091) q[3];
sx q[3];
rz(-1.3418875) q[3];
sx q[3];
rz(2.6218417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2733549) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-2.3256425) q[2];
rz(-0.50968918) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(1.9036487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8907392) q[0];
sx q[0];
rz(-1.7128523) q[0];
sx q[0];
rz(0.2302641) q[0];
rz(2.5157805) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(0.65840107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(-0.66396873) q[0];
x q[1];
rz(1.2344822) q[2];
sx q[2];
rz(-0.92421495) q[2];
sx q[2];
rz(0.99415776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2645996) q[1];
sx q[1];
rz(-2.0473285) q[1];
sx q[1];
rz(-0.18696733) q[1];
x q[2];
rz(-0.53604605) q[3];
sx q[3];
rz(-1.8744161) q[3];
sx q[3];
rz(-2.924078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-0.8224951) q[2];
sx q[2];
rz(-0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.7815536) q[3];
sx q[3];
rz(1.9256928) q[3];
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
rz(2.6170549) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(-1.9032003) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-0.17910437) q[2];
sx q[2];
rz(-2.1894107) q[2];
sx q[2];
rz(-1.4084963) q[2];
rz(1.5068549) q[3];
sx q[3];
rz(-0.99935525) q[3];
sx q[3];
rz(2.9730848) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
