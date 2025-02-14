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
rz(-1.8981847) q[0];
sx q[0];
rz(-1.564448) q[0];
sx q[0];
rz(2.2301883) q[0];
rz(2.4721594) q[1];
sx q[1];
rz(3.4179847) q[1];
sx q[1];
rz(8.5951947) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4938717) q[0];
sx q[0];
rz(-2.6135178) q[0];
sx q[0];
rz(-2.2654669) q[0];
rz(-0.24721036) q[2];
sx q[2];
rz(-1.6610378) q[2];
sx q[2];
rz(1.6177141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1948283) q[1];
sx q[1];
rz(-0.84673893) q[1];
sx q[1];
rz(1.9763234) q[1];
x q[2];
rz(-2.8220956) q[3];
sx q[3];
rz(-2.6852859) q[3];
sx q[3];
rz(2.1647477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7726045) q[2];
sx q[2];
rz(-0.96869865) q[2];
sx q[2];
rz(1.9770835) q[2];
rz(0.72689593) q[3];
sx q[3];
rz(-1.8098857) q[3];
sx q[3];
rz(-0.070778457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0182536) q[0];
sx q[0];
rz(-0.66272074) q[0];
sx q[0];
rz(0.29139274) q[0];
rz(0.60028752) q[1];
sx q[1];
rz(-2.184506) q[1];
sx q[1];
rz(0.16500638) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5614189) q[0];
sx q[0];
rz(-2.8695171) q[0];
sx q[0];
rz(1.838057) q[0];
x q[1];
rz(-1.8012456) q[2];
sx q[2];
rz(-2.682095) q[2];
sx q[2];
rz(0.27333096) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9143888) q[1];
sx q[1];
rz(-1.5014577) q[1];
sx q[1];
rz(0.68080618) q[1];
x q[2];
rz(0.21221186) q[3];
sx q[3];
rz(-2.8877875) q[3];
sx q[3];
rz(-0.55849822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0642285) q[2];
sx q[2];
rz(-0.95244971) q[2];
sx q[2];
rz(2.9020818) q[2];
rz(-2.8149878) q[3];
sx q[3];
rz(-2.9686847) q[3];
sx q[3];
rz(-3.0467564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8851682) q[0];
sx q[0];
rz(-0.435985) q[0];
sx q[0];
rz(-1.5727795) q[0];
rz(-0.39372152) q[1];
sx q[1];
rz(-2.5865793) q[1];
sx q[1];
rz(-0.47600019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9565292) q[0];
sx q[0];
rz(-1.0165043) q[0];
sx q[0];
rz(-2.8962062) q[0];
rz(-pi) q[1];
rz(-1.3838686) q[2];
sx q[2];
rz(-1.6003967) q[2];
sx q[2];
rz(1.6254978) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16533599) q[1];
sx q[1];
rz(-1.1603198) q[1];
sx q[1];
rz(2.2056376) q[1];
x q[2];
rz(0.94635238) q[3];
sx q[3];
rz(-0.40842429) q[3];
sx q[3];
rz(-2.8129336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3289045) q[2];
sx q[2];
rz(-0.71949553) q[2];
sx q[2];
rz(-2.1997814) q[2];
rz(-1.719126) q[3];
sx q[3];
rz(-2.7673281) q[3];
sx q[3];
rz(-2.4678738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0734237) q[0];
sx q[0];
rz(-0.24003679) q[0];
sx q[0];
rz(1.5330676) q[0];
rz(0.34670058) q[1];
sx q[1];
rz(-2.0997212) q[1];
sx q[1];
rz(-0.18992058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3118211) q[0];
sx q[0];
rz(-1.8320727) q[0];
sx q[0];
rz(-0.7423865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0937017) q[2];
sx q[2];
rz(-1.8162842) q[2];
sx q[2];
rz(2.7141822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0611113) q[1];
sx q[1];
rz(-0.60855344) q[1];
sx q[1];
rz(-1.4707277) q[1];
x q[2];
rz(-2.3307822) q[3];
sx q[3];
rz(-2.5805507) q[3];
sx q[3];
rz(-2.0138134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7803663) q[2];
sx q[2];
rz(-2.2915338) q[2];
sx q[2];
rz(-0.38193199) q[2];
rz(-3.0319038) q[3];
sx q[3];
rz(-2.1083125) q[3];
sx q[3];
rz(0.1194574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0383976) q[0];
sx q[0];
rz(-1.1829475) q[0];
sx q[0];
rz(0.84684816) q[0];
rz(3.0021744) q[1];
sx q[1];
rz(-2.3250695) q[1];
sx q[1];
rz(0.74904186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232305) q[0];
sx q[0];
rz(-1.4233755) q[0];
sx q[0];
rz(-0.69578232) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48495082) q[2];
sx q[2];
rz(-1.6963939) q[2];
sx q[2];
rz(1.6750977) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3917424) q[1];
sx q[1];
rz(-0.23255177) q[1];
sx q[1];
rz(-0.58156965) q[1];
rz(1.7827634) q[3];
sx q[3];
rz(-0.51111584) q[3];
sx q[3];
rz(-1.2455452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.9355115) q[2];
sx q[2];
rz(-1.2568018) q[2];
sx q[2];
rz(-1.067266) q[2];
rz(-2.4308128) q[3];
sx q[3];
rz(-1.658541) q[3];
sx q[3];
rz(1.425364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7988605) q[0];
sx q[0];
rz(-1.6628168) q[0];
sx q[0];
rz(-2.1730098) q[0];
rz(-2.4505278) q[1];
sx q[1];
rz(-1.7993118) q[1];
sx q[1];
rz(1.3074646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7846061) q[0];
sx q[0];
rz(-0.95080599) q[0];
sx q[0];
rz(1.8852573) q[0];
rz(-1.3771626) q[2];
sx q[2];
rz(-1.8566974) q[2];
sx q[2];
rz(1.4833956) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4279514) q[1];
sx q[1];
rz(-0.390807) q[1];
sx q[1];
rz(1.4321838) q[1];
rz(-pi) q[2];
rz(-1.8400165) q[3];
sx q[3];
rz(-1.3938187) q[3];
sx q[3];
rz(-0.38401803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4751733) q[2];
sx q[2];
rz(-0.23491493) q[2];
sx q[2];
rz(-1.2972181) q[2];
rz(-1.746486) q[3];
sx q[3];
rz(-1.3336983) q[3];
sx q[3];
rz(-1.6102128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.34201) q[0];
sx q[0];
rz(-0.81556773) q[0];
sx q[0];
rz(2.2014501) q[0];
rz(-0.6483342) q[1];
sx q[1];
rz(-1.2038566) q[1];
sx q[1];
rz(2.8727093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4057345) q[0];
sx q[0];
rz(-2.5395505) q[0];
sx q[0];
rz(-2.6989537) q[0];
rz(-1.2598632) q[2];
sx q[2];
rz(-2.249392) q[2];
sx q[2];
rz(2.1846365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2317048) q[1];
sx q[1];
rz(-2.1429166) q[1];
sx q[1];
rz(-0.27426274) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2327345) q[3];
sx q[3];
rz(-1.1629075) q[3];
sx q[3];
rz(3.0666587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1342643) q[2];
sx q[2];
rz(-1.2790044) q[2];
sx q[2];
rz(-2.4556665) q[2];
rz(0.736233) q[3];
sx q[3];
rz(-2.1015621) q[3];
sx q[3];
rz(-0.22549103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4418884) q[0];
sx q[0];
rz(-1.9206973) q[0];
sx q[0];
rz(0.7487444) q[0];
rz(0.092983149) q[1];
sx q[1];
rz(-2.3817606) q[1];
sx q[1];
rz(2.0704796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1919928) q[0];
sx q[0];
rz(-1.1566391) q[0];
sx q[0];
rz(-0.55121514) q[0];
rz(-pi) q[1];
rz(-2.4872141) q[2];
sx q[2];
rz(-1.1977695) q[2];
sx q[2];
rz(-1.4747064) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86102137) q[1];
sx q[1];
rz(-1.8771267) q[1];
sx q[1];
rz(2.2561399) q[1];
rz(-pi) q[2];
rz(-1.1597667) q[3];
sx q[3];
rz(-1.5229406) q[3];
sx q[3];
rz(1.4519297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45212713) q[2];
sx q[2];
rz(-2.5048544) q[2];
sx q[2];
rz(2.9316736) q[2];
rz(-0.12957761) q[3];
sx q[3];
rz(-2.1942873) q[3];
sx q[3];
rz(-1.5642222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.121948) q[0];
sx q[0];
rz(-0.26109281) q[0];
sx q[0];
rz(-0.013539465) q[0];
rz(2.6059222) q[1];
sx q[1];
rz(-0.56709138) q[1];
sx q[1];
rz(-0.013669107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5807728) q[0];
sx q[0];
rz(-2.3285064) q[0];
sx q[0];
rz(1.8845425) q[0];
x q[1];
rz(-2.7510468) q[2];
sx q[2];
rz(-1.1154113) q[2];
sx q[2];
rz(2.9519338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7175258) q[1];
sx q[1];
rz(-2.2560511) q[1];
sx q[1];
rz(-1.9937033) q[1];
rz(-1.9433653) q[3];
sx q[3];
rz(-2.2856224) q[3];
sx q[3];
rz(-1.4095588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6232264) q[2];
sx q[2];
rz(-1.7609111) q[2];
sx q[2];
rz(1.7784485) q[2];
rz(-2.3220883) q[3];
sx q[3];
rz(-0.47911152) q[3];
sx q[3];
rz(2.4578186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0226456) q[0];
sx q[0];
rz(-2.2830257) q[0];
sx q[0];
rz(-1.488142) q[0];
rz(0.34291521) q[1];
sx q[1];
rz(-1.5577134) q[1];
sx q[1];
rz(2.3905579) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2011178) q[0];
sx q[0];
rz(-1.6342666) q[0];
sx q[0];
rz(-1.8497784) q[0];
x q[1];
rz(-0.20047156) q[2];
sx q[2];
rz(-1.1159917) q[2];
sx q[2];
rz(-1.7998296) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4655059) q[1];
sx q[1];
rz(-1.6889531) q[1];
sx q[1];
rz(2.5398269) q[1];
rz(-1.6609825) q[3];
sx q[3];
rz(-1.5377561) q[3];
sx q[3];
rz(-0.15109135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.47120961) q[2];
sx q[2];
rz(-0.94380108) q[2];
sx q[2];
rz(2.75441) q[2];
rz(-0.47532982) q[3];
sx q[3];
rz(-1.9742222) q[3];
sx q[3];
rz(2.3502684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925856) q[0];
sx q[0];
rz(-2.1795166) q[0];
sx q[0];
rz(-2.4695061) q[0];
rz(2.7723906) q[1];
sx q[1];
rz(-2.3875356) q[1];
sx q[1];
rz(1.748132) q[1];
rz(2.4245928) q[2];
sx q[2];
rz(-1.430027) q[2];
sx q[2];
rz(1.2363557) q[2];
rz(1.297507) q[3];
sx q[3];
rz(-1.6798301) q[3];
sx q[3];
rz(-1.9435431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
