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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(-2.1313957) q[0];
rz(-2.9695192) q[1];
sx q[1];
rz(-0.80614027) q[1];
sx q[1];
rz(-2.1176718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14903325) q[0];
sx q[0];
rz(-1.5776538) q[0];
sx q[0];
rz(-2.0586781) q[0];
x q[1];
rz(2.6226343) q[2];
sx q[2];
rz(-0.4441688) q[2];
sx q[2];
rz(-0.78625317) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5883903) q[1];
sx q[1];
rz(-1.4103048) q[1];
sx q[1];
rz(-3.1218916) q[1];
rz(1.6590622) q[3];
sx q[3];
rz(-1.3219943) q[3];
sx q[3];
rz(-0.036347957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94850928) q[2];
sx q[2];
rz(-2.812959) q[2];
sx q[2];
rz(2.9060717) q[2];
rz(2.113302) q[3];
sx q[3];
rz(-1.8833269) q[3];
sx q[3];
rz(-1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368206) q[0];
sx q[0];
rz(-1.7828159) q[0];
sx q[0];
rz(2.219668) q[0];
rz(3.0251265) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(-2.2540653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1257) q[0];
sx q[0];
rz(-1.757048) q[0];
sx q[0];
rz(-0.13293477) q[0];
x q[1];
rz(-2.791128) q[2];
sx q[2];
rz(-1.6258569) q[2];
sx q[2];
rz(0.45086917) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16819363) q[1];
sx q[1];
rz(-2.138294) q[1];
sx q[1];
rz(1.7048111) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87882249) q[3];
sx q[3];
rz(-1.678595) q[3];
sx q[3];
rz(-2.381058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5827215) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(0.23898807) q[2];
rz(-0.72238266) q[3];
sx q[3];
rz(-0.84010774) q[3];
sx q[3];
rz(1.9766138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8257985) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(2.353299) q[0];
rz(1.7514508) q[1];
sx q[1];
rz(-0.43594053) q[1];
sx q[1];
rz(1.4424666) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4460488) q[0];
sx q[0];
rz(-0.15803495) q[0];
sx q[0];
rz(-1.1084854) q[0];
rz(-pi) q[1];
x q[1];
rz(1.945247) q[2];
sx q[2];
rz(-1.2948639) q[2];
sx q[2];
rz(2.3464613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6953814) q[1];
sx q[1];
rz(-1.808481) q[1];
sx q[1];
rz(-0.85376815) q[1];
rz(-pi) q[2];
rz(1.0679108) q[3];
sx q[3];
rz(-1.4137795) q[3];
sx q[3];
rz(-2.5105421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.73551377) q[2];
sx q[2];
rz(-1.1512076) q[2];
sx q[2];
rz(2.8847412) q[2];
rz(-0.86366051) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(-2.5679722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470806) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(-1.5091913) q[0];
rz(-2.8250561) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(1.0155771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7596801) q[0];
sx q[0];
rz(-2.0354871) q[0];
sx q[0];
rz(-1.0601935) q[0];
rz(-pi) q[1];
rz(-0.79358856) q[2];
sx q[2];
rz(-1.3745135) q[2];
sx q[2];
rz(-1.0283141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16962651) q[1];
sx q[1];
rz(-1.7952654) q[1];
sx q[1];
rz(-1.9588934) q[1];
rz(-1.7523132) q[3];
sx q[3];
rz(-0.75541516) q[3];
sx q[3];
rz(0.28873539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(1.451937) q[2];
rz(1.9076094) q[3];
sx q[3];
rz(-2.0472287) q[3];
sx q[3];
rz(0.88172495) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7538309) q[0];
sx q[0];
rz(-1.3085288) q[0];
sx q[0];
rz(-1.0863139) q[0];
rz(-1.4007252) q[1];
sx q[1];
rz(-2.3298405) q[1];
sx q[1];
rz(-3.1382255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683748) q[0];
sx q[0];
rz(-1.5766607) q[0];
sx q[0];
rz(0.59451367) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1589627) q[2];
sx q[2];
rz(-1.8525043) q[2];
sx q[2];
rz(-2.1717193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7351384) q[1];
sx q[1];
rz(-0.5105831) q[1];
sx q[1];
rz(-0.11514727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5967054) q[3];
sx q[3];
rz(-2.1707275) q[3];
sx q[3];
rz(-2.5480218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1668261) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(-0.33494803) q[2];
rz(-2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(-3.0193442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5122546) q[0];
sx q[0];
rz(-0.4991931) q[0];
sx q[0];
rz(2.7235624) q[0];
rz(2.4755075) q[1];
sx q[1];
rz(-1.8981551) q[1];
sx q[1];
rz(-2.8935249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6757096) q[0];
sx q[0];
rz(-1.0843236) q[0];
sx q[0];
rz(-2.7652362) q[0];
x q[1];
rz(1.4573967) q[2];
sx q[2];
rz(-2.6893986) q[2];
sx q[2];
rz(2.3101286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0736165) q[1];
sx q[1];
rz(-2.4779137) q[1];
sx q[1];
rz(-0.90866088) q[1];
rz(-pi) q[2];
rz(-0.25603981) q[3];
sx q[3];
rz(-2.1845792) q[3];
sx q[3];
rz(-0.41071329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4416113) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(-1.1629533) q[2];
rz(-0.59398389) q[3];
sx q[3];
rz(-1.5409527) q[3];
sx q[3];
rz(-0.99015132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6043337) q[0];
sx q[0];
rz(-0.27892497) q[0];
sx q[0];
rz(3.0766686) q[0];
rz(0.36554947) q[1];
sx q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(-2.4117267) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.192963) q[0];
sx q[0];
rz(-2.2092487) q[0];
sx q[0];
rz(3.0045511) q[0];
rz(-pi) q[1];
rz(1.7634298) q[2];
sx q[2];
rz(-2.3245077) q[2];
sx q[2];
rz(1.1258923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0522451) q[1];
sx q[1];
rz(-2.5762853) q[1];
sx q[1];
rz(1.7491399) q[1];
rz(-pi) q[2];
rz(2.6486126) q[3];
sx q[3];
rz(-2.3928309) q[3];
sx q[3];
rz(-1.6822588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2713752) q[2];
sx q[2];
rz(-1.368618) q[2];
sx q[2];
rz(-1.082487) q[2];
rz(-1.4879976) q[3];
sx q[3];
rz(-2.7797785) q[3];
sx q[3];
rz(-2.7910119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(0.3311232) q[0];
rz(1.6638727) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(0.32041916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54658157) q[0];
sx q[0];
rz(-3.0394331) q[0];
sx q[0];
rz(2.8861396) q[0];
x q[1];
rz(-1.7170948) q[2];
sx q[2];
rz(-1.4182036) q[2];
sx q[2];
rz(-2.1426147) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72057331) q[1];
sx q[1];
rz(-1.5552727) q[1];
sx q[1];
rz(0.42625721) q[1];
rz(-1.5005169) q[3];
sx q[3];
rz(-2.456074) q[3];
sx q[3];
rz(0.34834114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42935818) q[2];
sx q[2];
rz(-1.7244491) q[2];
sx q[2];
rz(-2.7799535) q[2];
rz(-2.2278191) q[3];
sx q[3];
rz(-0.29699609) q[3];
sx q[3];
rz(0.44338068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2509895) q[0];
sx q[0];
rz(-2.3539703) q[0];
sx q[0];
rz(-3.0255764) q[0];
rz(-1.8138255) q[1];
sx q[1];
rz(-1.9783744) q[1];
sx q[1];
rz(1.2001002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2105149) q[0];
sx q[0];
rz(-1.7621982) q[0];
sx q[0];
rz(1.5516419) q[0];
x q[1];
rz(-2.2564327) q[2];
sx q[2];
rz(-1.6994393) q[2];
sx q[2];
rz(2.7264915) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45428206) q[1];
sx q[1];
rz(-0.56182278) q[1];
sx q[1];
rz(0.027472036) q[1];
rz(2.9755244) q[3];
sx q[3];
rz(-1.3343873) q[3];
sx q[3];
rz(-2.0061559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3001331) q[2];
sx q[2];
rz(-1.4347142) q[2];
sx q[2];
rz(2.8099698) q[2];
rz(-0.2837818) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(0.42154977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3109741) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(-2.9750138) q[0];
rz(2.2919948) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(3.0679682) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1243116) q[0];
sx q[0];
rz(-0.86605427) q[0];
sx q[0];
rz(1.2344735) q[0];
x q[1];
rz(2.2715015) q[2];
sx q[2];
rz(-2.4656762) q[2];
sx q[2];
rz(1.5004339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6839809) q[1];
sx q[1];
rz(-1.906412) q[1];
sx q[1];
rz(2.6102351) q[1];
rz(-1.0898548) q[3];
sx q[3];
rz(-2.1409799) q[3];
sx q[3];
rz(2.7358304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9061884) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(-0.6582312) q[2];
rz(1.8002347) q[3];
sx q[3];
rz(-0.96786371) q[3];
sx q[3];
rz(-2.2906176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8258719) q[0];
sx q[0];
rz(-2.3128339) q[0];
sx q[0];
rz(-0.88902938) q[0];
rz(-0.52147621) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(0.78962959) q[2];
sx q[2];
rz(-1.9676002) q[2];
sx q[2];
rz(-1.0330539) q[2];
rz(-1.7375199) q[3];
sx q[3];
rz(-2.1862595) q[3];
sx q[3];
rz(-1.4487472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
