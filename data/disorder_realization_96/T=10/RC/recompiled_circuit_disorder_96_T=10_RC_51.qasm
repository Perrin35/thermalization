OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6450206) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(-0.42005959) q[0];
rz(-pi) q[1];
rz(2.1098233) q[2];
sx q[2];
rz(-1.6845778) q[2];
sx q[2];
rz(0.10345085) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9125036) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(1.9181262) q[1];
x q[2];
rz(0.79521631) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(0.29602805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(-1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5553404) q[0];
sx q[0];
rz(-2.147701) q[0];
sx q[0];
rz(0.14848407) q[0];
rz(1.8230121) q[2];
sx q[2];
rz(-2.2644342) q[2];
sx q[2];
rz(-1.1064305) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10837308) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(0.56652041) q[1];
x q[2];
rz(-1.4278533) q[3];
sx q[3];
rz(-1.6030451) q[3];
sx q[3];
rz(-1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0216996) q[2];
sx q[2];
rz(-1.0571486) q[2];
sx q[2];
rz(-0.95834857) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(-2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2496101) q[0];
sx q[0];
rz(-1.550807) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
x q[1];
rz(1.259272) q[2];
sx q[2];
rz(-1.4853962) q[2];
sx q[2];
rz(0.66750079) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5379996) q[1];
sx q[1];
rz(-2.2003761) q[1];
sx q[1];
rz(-1.0521207) q[1];
x q[2];
rz(-1.2529536) q[3];
sx q[3];
rz(-0.40099537) q[3];
sx q[3];
rz(-2.3272115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(0.5853816) q[2];
rz(0.18150005) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25585184) q[0];
sx q[0];
rz(-0.8937853) q[0];
sx q[0];
rz(2.6141502) q[0];
rz(-2.9085607) q[2];
sx q[2];
rz(-0.2430025) q[2];
sx q[2];
rz(2.255893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.199898) q[1];
sx q[1];
rz(-2.3823793) q[1];
sx q[1];
rz(0.56337507) q[1];
x q[2];
rz(-0.67646497) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(1.3456618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(2.0902324) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(0.043118127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2152104) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-0.2290639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28508614) q[2];
sx q[2];
rz(-1.0694155) q[2];
sx q[2];
rz(1.8269055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4700714) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(2.3684711) q[1];
x q[2];
rz(-0.78955663) q[3];
sx q[3];
rz(-1.4649179) q[3];
sx q[3];
rz(1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0126426) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(2.5514065) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689894) q[0];
sx q[0];
rz(-0.4840584) q[0];
sx q[0];
rz(2.6812535) q[0];
x q[1];
rz(1.9361587) q[2];
sx q[2];
rz(-2.0632671) q[2];
sx q[2];
rz(1.9838711) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(0.92909716) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70763564) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(-1.9572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6713312) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0253042) q[0];
sx q[0];
rz(-1.7968654) q[0];
sx q[0];
rz(-2.6095005) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5146671) q[2];
sx q[2];
rz(-1.727384) q[2];
sx q[2];
rz(2.2121034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0824273) q[1];
sx q[1];
rz(-1.6470243) q[1];
sx q[1];
rz(0.04870292) q[1];
rz(-pi) q[2];
rz(2.9229786) q[3];
sx q[3];
rz(-0.84158763) q[3];
sx q[3];
rz(1.5515755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-0.41762525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324875) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(-2.0672654) q[0];
rz(-pi) q[1];
rz(-3.0268961) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(2.7483658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3216746) q[1];
sx q[1];
rz(-2.3707317) q[1];
sx q[1];
rz(2.2714771) q[1];
rz(-pi) q[2];
rz(-0.91248625) q[3];
sx q[3];
rz(-0.45966002) q[3];
sx q[3];
rz(1.0115136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8326571) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(0.94690698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.2456018) q[0];
sx q[0];
rz(1.253771) q[0];
rz(1.6106748) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(-3.1181042) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.088989181) q[1];
sx q[1];
rz(-2.7967884) q[1];
sx q[1];
rz(0.25119541) q[1];
x q[2];
rz(-2.0128653) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(-2.5069619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56069121) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(2.4272264) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(-2.4972829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2243758) q[0];
sx q[0];
rz(-1.7443568) q[0];
sx q[0];
rz(1.5149084) q[0];
rz(-2.2868025) q[2];
sx q[2];
rz(-1.3959179) q[2];
sx q[2];
rz(-2.1388432) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.555968) q[1];
sx q[1];
rz(-1.854419) q[1];
sx q[1];
rz(2.7908299) q[1];
rz(-pi) q[2];
rz(-2.042291) q[3];
sx q[3];
rz(-1.9715371) q[3];
sx q[3];
rz(2.5221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(-2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-1.3953801) q[2];
sx q[2];
rz(-2.5098364) q[2];
sx q[2];
rz(-2.9927158) q[2];
rz(-2.2453528) q[3];
sx q[3];
rz(-1.2590209) q[3];
sx q[3];
rz(-0.38929064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
