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
rz(-1.9528376) q[0];
sx q[0];
rz(-0.83851695) q[0];
sx q[0];
rz(-1.4128348) q[0];
rz(-2.7589594) q[1];
sx q[1];
rz(-1.1525947) q[1];
sx q[1];
rz(-1.1421854) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4948311) q[0];
sx q[0];
rz(-1.9660304) q[0];
sx q[0];
rz(-3.0656205) q[0];
x q[1];
rz(0.5025592) q[2];
sx q[2];
rz(-1.1389073) q[2];
sx q[2];
rz(2.2186861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22622977) q[1];
sx q[1];
rz(-2.1945476) q[1];
sx q[1];
rz(2.8585494) q[1];
x q[2];
rz(0.19742404) q[3];
sx q[3];
rz(-0.89251489) q[3];
sx q[3];
rz(-1.9367644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4932058) q[2];
sx q[2];
rz(-1.9923261) q[2];
sx q[2];
rz(-3.1374078) q[2];
rz(-0.74797136) q[3];
sx q[3];
rz(-2.3145521) q[3];
sx q[3];
rz(0.23676693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10906049) q[0];
sx q[0];
rz(-2.2219658) q[0];
sx q[0];
rz(0.064706651) q[0];
rz(1.0626571) q[1];
sx q[1];
rz(-0.92667842) q[1];
sx q[1];
rz(1.0437171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1239584) q[0];
sx q[0];
rz(-1.7746266) q[0];
sx q[0];
rz(-2.6685326) q[0];
x q[1];
rz(-0.02702464) q[2];
sx q[2];
rz(-1.6690429) q[2];
sx q[2];
rz(2.8842928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1295075) q[1];
sx q[1];
rz(-1.7237067) q[1];
sx q[1];
rz(2.6430886) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6297335) q[3];
sx q[3];
rz(-2.8843237) q[3];
sx q[3];
rz(2.9016837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8740497) q[2];
sx q[2];
rz(-1.1288613) q[2];
sx q[2];
rz(-1.4675325) q[2];
rz(2.2360146) q[3];
sx q[3];
rz(-1.0970683) q[3];
sx q[3];
rz(0.21152285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.31731376) q[0];
sx q[0];
rz(-0.86238328) q[0];
sx q[0];
rz(2.8381919) q[0];
rz(0.41269451) q[1];
sx q[1];
rz(-1.7957325) q[1];
sx q[1];
rz(-2.5491098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68776206) q[0];
sx q[0];
rz(-2.651281) q[0];
sx q[0];
rz(-0.47435905) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7492152) q[2];
sx q[2];
rz(-2.1270299) q[2];
sx q[2];
rz(-1.2715593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35960782) q[1];
sx q[1];
rz(-2.4949412) q[1];
sx q[1];
rz(-2.9293865) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6487916) q[3];
sx q[3];
rz(-0.78712026) q[3];
sx q[3];
rz(-1.9446179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5182284) q[2];
sx q[2];
rz(-2.6609504) q[2];
sx q[2];
rz(-0.19409689) q[2];
rz(0.57191166) q[3];
sx q[3];
rz(-1.5587991) q[3];
sx q[3];
rz(1.5863878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4647575) q[0];
sx q[0];
rz(-1.06523) q[0];
sx q[0];
rz(1.42365) q[0];
rz(-0.33547297) q[1];
sx q[1];
rz(-2.5392541) q[1];
sx q[1];
rz(2.4981892) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5945054) q[0];
sx q[0];
rz(-2.3427561) q[0];
sx q[0];
rz(-0.79669768) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5667666) q[2];
sx q[2];
rz(-1.4919623) q[2];
sx q[2];
rz(-3.0835033) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5550892) q[1];
sx q[1];
rz(-2.2528306) q[1];
sx q[1];
rz(-2.5413496) q[1];
rz(-1.8845064) q[3];
sx q[3];
rz(-2.6711552) q[3];
sx q[3];
rz(1.1251118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.246835) q[2];
sx q[2];
rz(-1.094123) q[2];
sx q[2];
rz(0.82478729) q[2];
rz(0.60683933) q[3];
sx q[3];
rz(-1.9347128) q[3];
sx q[3];
rz(-1.3185893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987228) q[0];
sx q[0];
rz(-0.56951183) q[0];
sx q[0];
rz(-0.2670162) q[0];
rz(0.46785242) q[1];
sx q[1];
rz(-2.0303969) q[1];
sx q[1];
rz(1.1763447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699791) q[0];
sx q[0];
rz(-1.5494657) q[0];
sx q[0];
rz(1.5202549) q[0];
rz(-0.77645923) q[2];
sx q[2];
rz(-1.0506223) q[2];
sx q[2];
rz(0.3683683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5084978) q[1];
sx q[1];
rz(-1.1207241) q[1];
sx q[1];
rz(-3.0799559) q[1];
x q[2];
rz(2.4202706) q[3];
sx q[3];
rz(-1.112794) q[3];
sx q[3];
rz(1.3617226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.10845575) q[2];
sx q[2];
rz(-1.5908073) q[2];
sx q[2];
rz(0.94672686) q[2];
rz(1.4416134) q[3];
sx q[3];
rz(-1.0333034) q[3];
sx q[3];
rz(1.7472964) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93447584) q[0];
sx q[0];
rz(-1.8659135) q[0];
sx q[0];
rz(1.3502655) q[0];
rz(2.4840202) q[1];
sx q[1];
rz(-1.5788199) q[1];
sx q[1];
rz(-1.879976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3917106) q[0];
sx q[0];
rz(-0.80750033) q[0];
sx q[0];
rz(0.9901643) q[0];
x q[1];
rz(0.96716934) q[2];
sx q[2];
rz(-2.1534922) q[2];
sx q[2];
rz(0.36925579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2985136) q[1];
sx q[1];
rz(-1.635478) q[1];
sx q[1];
rz(-2.9119125) q[1];
rz(1.3587977) q[3];
sx q[3];
rz(-1.7431211) q[3];
sx q[3];
rz(1.621877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6235846) q[2];
sx q[2];
rz(-1.1058747) q[2];
sx q[2];
rz(0.27274954) q[2];
rz(-2.2128211) q[3];
sx q[3];
rz(-0.85799587) q[3];
sx q[3];
rz(1.4906918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51468325) q[0];
sx q[0];
rz(-2.1105483) q[0];
sx q[0];
rz(-0.92460257) q[0];
rz(0.54479105) q[1];
sx q[1];
rz(-1.7926615) q[1];
sx q[1];
rz(-3.0371688) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6729148) q[0];
sx q[0];
rz(-0.12355655) q[0];
sx q[0];
rz(-2.3455023) q[0];
x q[1];
rz(-2.7265276) q[2];
sx q[2];
rz(-1.7156148) q[2];
sx q[2];
rz(-0.96127779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.50867289) q[1];
sx q[1];
rz(-1.3471222) q[1];
sx q[1];
rz(-1.2583076) q[1];
rz(-pi) q[2];
rz(-0.35108836) q[3];
sx q[3];
rz(-1.2005873) q[3];
sx q[3];
rz(2.5552487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8203848) q[2];
sx q[2];
rz(-1.2106004) q[2];
sx q[2];
rz(-1.661181) q[2];
rz(3.1373451) q[3];
sx q[3];
rz(-0.85631266) q[3];
sx q[3];
rz(2.492304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5874807) q[0];
sx q[0];
rz(-0.09621796) q[0];
sx q[0];
rz(-1.2531248) q[0];
rz(-2.9907277) q[1];
sx q[1];
rz(-0.8395218) q[1];
sx q[1];
rz(1.8419267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74053451) q[0];
sx q[0];
rz(-0.63176934) q[0];
sx q[0];
rz(-0.44332544) q[0];
x q[1];
rz(1.9076882) q[2];
sx q[2];
rz(-0.46122069) q[2];
sx q[2];
rz(2.5835844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7902657) q[1];
sx q[1];
rz(-0.68584397) q[1];
sx q[1];
rz(-3.1180361) q[1];
rz(2.9392936) q[3];
sx q[3];
rz(-1.222156) q[3];
sx q[3];
rz(-1.543247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3018939) q[2];
sx q[2];
rz(-1.7599186) q[2];
sx q[2];
rz(1.1559486) q[2];
rz(2.4902952) q[3];
sx q[3];
rz(-0.39172253) q[3];
sx q[3];
rz(2.7294066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10525178) q[0];
sx q[0];
rz(-1.3407433) q[0];
sx q[0];
rz(0.80129188) q[0];
rz(0.85805145) q[1];
sx q[1];
rz(-0.67887226) q[1];
sx q[1];
rz(-0.37989315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44110041) q[0];
sx q[0];
rz(-2.0394562) q[0];
sx q[0];
rz(2.9392713) q[0];
x q[1];
rz(1.395389) q[2];
sx q[2];
rz(-1.6885969) q[2];
sx q[2];
rz(0.30994383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2945248) q[1];
sx q[1];
rz(-2.2502568) q[1];
sx q[1];
rz(-0.39107283) q[1];
rz(-3.1201023) q[3];
sx q[3];
rz(-2.700319) q[3];
sx q[3];
rz(3.0336371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51565591) q[2];
sx q[2];
rz(-0.3346602) q[2];
sx q[2];
rz(-2.3626309) q[2];
rz(-0.37911478) q[3];
sx q[3];
rz(-1.4232676) q[3];
sx q[3];
rz(0.085722119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871386) q[0];
sx q[0];
rz(-1.4105281) q[0];
sx q[0];
rz(-2.6692303) q[0];
rz(-2.8624599) q[1];
sx q[1];
rz(-2.0886853) q[1];
sx q[1];
rz(0.36769029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5147509) q[0];
sx q[0];
rz(-1.923133) q[0];
sx q[0];
rz(0.3120089) q[0];
x q[1];
rz(-0.49007551) q[2];
sx q[2];
rz(-0.44675175) q[2];
sx q[2];
rz(1.2778145) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.20357) q[1];
sx q[1];
rz(-1.1633151) q[1];
sx q[1];
rz(-1.8304871) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84334737) q[3];
sx q[3];
rz(-0.9002004) q[3];
sx q[3];
rz(1.5274109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3231861) q[2];
sx q[2];
rz(-1.0419934) q[2];
sx q[2];
rz(-0.47427487) q[2];
rz(0.91931528) q[3];
sx q[3];
rz(-1.5761458) q[3];
sx q[3];
rz(-1.7491755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11488386) q[0];
sx q[0];
rz(-1.8395431) q[0];
sx q[0];
rz(-1.607847) q[0];
rz(-2.907091) q[1];
sx q[1];
rz(-1.0954183) q[1];
sx q[1];
rz(-0.29874994) q[1];
rz(-0.99954323) q[2];
sx q[2];
rz(-1.4951946) q[2];
sx q[2];
rz(-1.9517938) q[2];
rz(1.5212223) q[3];
sx q[3];
rz(-1.5885316) q[3];
sx q[3];
rz(1.9902609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
