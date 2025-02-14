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
rz(-1.0722941) q[0];
sx q[0];
rz(-0.18360734) q[0];
sx q[0];
rz(-0.16820964) q[0];
rz(-1.3104562) q[1];
sx q[1];
rz(-2.0252731) q[1];
sx q[1];
rz(2.0235553) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71336761) q[0];
sx q[0];
rz(-1.8963739) q[0];
sx q[0];
rz(-1.6859562) q[0];
rz(-pi) q[1];
rz(-0.69324915) q[2];
sx q[2];
rz(-1.8680671) q[2];
sx q[2];
rz(1.3849486) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1217482) q[1];
sx q[1];
rz(-0.81070527) q[1];
sx q[1];
rz(1.5365063) q[1];
x q[2];
rz(-1.2994475) q[3];
sx q[3];
rz(-1.4381471) q[3];
sx q[3];
rz(2.1675183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85718021) q[2];
sx q[2];
rz(-0.99049157) q[2];
sx q[2];
rz(1.277479) q[2];
rz(-2.5340581) q[3];
sx q[3];
rz(-1.3527063) q[3];
sx q[3];
rz(-1.386462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69433895) q[0];
sx q[0];
rz(-0.67830938) q[0];
sx q[0];
rz(0.69764486) q[0];
rz(-2.8088226) q[1];
sx q[1];
rz(-2.0336626) q[1];
sx q[1];
rz(2.8093801) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32658122) q[0];
sx q[0];
rz(-1.7012811) q[0];
sx q[0];
rz(-1.4490118) q[0];
x q[1];
rz(-1.2142014) q[2];
sx q[2];
rz(-0.60672543) q[2];
sx q[2];
rz(2.1347097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4302286) q[1];
sx q[1];
rz(-1.9875184) q[1];
sx q[1];
rz(0.27208379) q[1];
rz(-2.6744889) q[3];
sx q[3];
rz(-2.3574266) q[3];
sx q[3];
rz(2.2968366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.701) q[2];
sx q[2];
rz(-1.168246) q[2];
sx q[2];
rz(-2.8458703) q[2];
rz(0.010585636) q[3];
sx q[3];
rz(-0.9897832) q[3];
sx q[3];
rz(-0.73339614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4877141) q[0];
sx q[0];
rz(-0.45844498) q[0];
sx q[0];
rz(2.606875) q[0];
rz(2.0405105) q[1];
sx q[1];
rz(-1.4357166) q[1];
sx q[1];
rz(-0.91317493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064414211) q[0];
sx q[0];
rz(-2.065669) q[0];
sx q[0];
rz(2.7735308) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9703254) q[2];
sx q[2];
rz(-2.5472884) q[2];
sx q[2];
rz(0.28501302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6102741) q[1];
sx q[1];
rz(-2.2372359) q[1];
sx q[1];
rz(-0.6998723) q[1];
rz(-0.072973386) q[3];
sx q[3];
rz(-2.2984031) q[3];
sx q[3];
rz(0.062372717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5494626) q[2];
sx q[2];
rz(-0.89207804) q[2];
sx q[2];
rz(-2.7590052) q[2];
rz(1.654918) q[3];
sx q[3];
rz(-2.8463709) q[3];
sx q[3];
rz(-2.2216643) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53052467) q[0];
sx q[0];
rz(-1.684573) q[0];
sx q[0];
rz(-1.0462421) q[0];
rz(-2.5776082) q[1];
sx q[1];
rz(-1.0020741) q[1];
sx q[1];
rz(1.3169588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038669) q[0];
sx q[0];
rz(-0.90529862) q[0];
sx q[0];
rz(-2.4212877) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1649873) q[2];
sx q[2];
rz(-0.69402678) q[2];
sx q[2];
rz(-2.9846689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0610132) q[1];
sx q[1];
rz(-1.6135585) q[1];
sx q[1];
rz(-0.31400915) q[1];
rz(-pi) q[2];
rz(-0.78329194) q[3];
sx q[3];
rz(-1.8785333) q[3];
sx q[3];
rz(0.35959343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4468533) q[2];
sx q[2];
rz(-1.6971735) q[2];
sx q[2];
rz(-2.1430338) q[2];
rz(0.57991943) q[3];
sx q[3];
rz(-1.4812352) q[3];
sx q[3];
rz(1.3045605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.21438433) q[0];
sx q[0];
rz(-1.2283607) q[0];
sx q[0];
rz(0.71980113) q[0];
rz(0.7431227) q[1];
sx q[1];
rz(-1.9235976) q[1];
sx q[1];
rz(-0.057083759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2190773) q[0];
sx q[0];
rz(-1.9382028) q[0];
sx q[0];
rz(1.7149345) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.916831) q[2];
sx q[2];
rz(-0.84362307) q[2];
sx q[2];
rz(2.7131701) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.067043153) q[1];
sx q[1];
rz(-0.32575575) q[1];
sx q[1];
rz(0.52149741) q[1];
x q[2];
rz(-2.1716398) q[3];
sx q[3];
rz(-0.68884736) q[3];
sx q[3];
rz(-0.057202489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7032949) q[2];
sx q[2];
rz(-1.1763828) q[2];
sx q[2];
rz(-1.8717742) q[2];
rz(0.54369175) q[3];
sx q[3];
rz(-0.68259493) q[3];
sx q[3];
rz(-1.3683176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0755587) q[0];
sx q[0];
rz(-2.7664001) q[0];
sx q[0];
rz(2.913108) q[0];
rz(0.98555073) q[1];
sx q[1];
rz(-1.5864213) q[1];
sx q[1];
rz(1.2062581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0096480308) q[0];
sx q[0];
rz(-1.7382657) q[0];
sx q[0];
rz(0.33958774) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2007063) q[2];
sx q[2];
rz(-0.36782757) q[2];
sx q[2];
rz(1.5076249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9555743) q[1];
sx q[1];
rz(-1.4494872) q[1];
sx q[1];
rz(-0.96069161) q[1];
rz(-pi) q[2];
rz(0.1749707) q[3];
sx q[3];
rz(-2.5911883) q[3];
sx q[3];
rz(-0.016684858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1365405) q[2];
sx q[2];
rz(-1.8661934) q[2];
sx q[2];
rz(1.8955815) q[2];
rz(-0.15240845) q[3];
sx q[3];
rz(-2.992575) q[3];
sx q[3];
rz(-2.3685031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8361255) q[0];
sx q[0];
rz(-1.9274599) q[0];
sx q[0];
rz(0.25492302) q[0];
rz(2.1615324) q[1];
sx q[1];
rz(-1.7995116) q[1];
sx q[1];
rz(0.23475383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7388044) q[0];
sx q[0];
rz(-2.1128034) q[0];
sx q[0];
rz(2.8554585) q[0];
rz(-pi) q[1];
rz(0.94257109) q[2];
sx q[2];
rz(-2.8543575) q[2];
sx q[2];
rz(0.38610215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3666881) q[1];
sx q[1];
rz(-1.8340381) q[1];
sx q[1];
rz(2.7145492) q[1];
x q[2];
rz(3.023769) q[3];
sx q[3];
rz(-2.243968) q[3];
sx q[3];
rz(-2.4073659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0945101) q[2];
sx q[2];
rz(-2.7213056) q[2];
sx q[2];
rz(-1.1959929) q[2];
rz(0.99810537) q[3];
sx q[3];
rz(-1.2752897) q[3];
sx q[3];
rz(2.3534145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274479) q[0];
sx q[0];
rz(-0.87004167) q[0];
sx q[0];
rz(-2.5092464) q[0];
rz(-1.5517722) q[1];
sx q[1];
rz(-1.3464144) q[1];
sx q[1];
rz(1.2987035) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.08399) q[0];
sx q[0];
rz(-2.6267536) q[0];
sx q[0];
rz(-1.043651) q[0];
rz(-0.66565158) q[2];
sx q[2];
rz(-1.5558262) q[2];
sx q[2];
rz(-0.41757181) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3159193) q[1];
sx q[1];
rz(-0.50858595) q[1];
sx q[1];
rz(3.0561165) q[1];
rz(-pi) q[2];
rz(-1.5428144) q[3];
sx q[3];
rz(-1.4194427) q[3];
sx q[3];
rz(1.7773113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87045264) q[2];
sx q[2];
rz(-2.2282034) q[2];
sx q[2];
rz(2.3823628) q[2];
rz(2.7325654) q[3];
sx q[3];
rz(-1.4141915) q[3];
sx q[3];
rz(3.0294688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2893534) q[0];
sx q[0];
rz(-1.2989346) q[0];
sx q[0];
rz(0.36561832) q[0];
rz(3.1386555) q[1];
sx q[1];
rz(-1.6117088) q[1];
sx q[1];
rz(-1.0681577) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.070052) q[0];
sx q[0];
rz(-2.5134633) q[0];
sx q[0];
rz(2.6578636) q[0];
rz(0.67068213) q[2];
sx q[2];
rz(-0.95206958) q[2];
sx q[2];
rz(2.3681792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5369027) q[1];
sx q[1];
rz(-0.3893756) q[1];
sx q[1];
rz(1.6714208) q[1];
x q[2];
rz(-0.2016507) q[3];
sx q[3];
rz(-1.5439873) q[3];
sx q[3];
rz(0.22814116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.7630161) q[2];
sx q[2];
rz(-0.91332674) q[2];
sx q[2];
rz(0.98706377) q[2];
rz(2.3297564) q[3];
sx q[3];
rz(-2.2631009) q[3];
sx q[3];
rz(1.2442376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064747485) q[0];
sx q[0];
rz(-0.61926121) q[0];
sx q[0];
rz(-1.6774119) q[0];
rz(2.1773416) q[1];
sx q[1];
rz(-1.3136656) q[1];
sx q[1];
rz(-1.3826694) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89870533) q[0];
sx q[0];
rz(-0.90685493) q[0];
sx q[0];
rz(2.6580145) q[0];
x q[1];
rz(-2.7708762) q[2];
sx q[2];
rz(-1.3138736) q[2];
sx q[2];
rz(-0.020250016) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9082038) q[1];
sx q[1];
rz(-2.1022132) q[1];
sx q[1];
rz(-2.6446093) q[1];
rz(-pi) q[2];
rz(-0.55158864) q[3];
sx q[3];
rz(-1.0304215) q[3];
sx q[3];
rz(2.0575749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1796403) q[2];
sx q[2];
rz(-1.3397168) q[2];
sx q[2];
rz(-1.4947653) q[2];
rz(-2.0252114) q[3];
sx q[3];
rz(-2.3514082) q[3];
sx q[3];
rz(-0.60379973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2434926) q[0];
sx q[0];
rz(-1.6742764) q[0];
sx q[0];
rz(1.8754638) q[0];
rz(0.19337868) q[1];
sx q[1];
rz(-2.0722176) q[1];
sx q[1];
rz(-1.8373012) q[1];
rz(1.8776352) q[2];
sx q[2];
rz(-2.2813792) q[2];
sx q[2];
rz(1.1165237) q[2];
rz(0.51283097) q[3];
sx q[3];
rz(-0.75779946) q[3];
sx q[3];
rz(2.1661314) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
