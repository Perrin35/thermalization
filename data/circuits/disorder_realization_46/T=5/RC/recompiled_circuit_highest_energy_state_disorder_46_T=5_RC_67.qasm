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
rz(1.8019851) q[0];
sx q[0];
rz(-2.7924502) q[0];
sx q[0];
rz(-2.1139297) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(-1.0171913) q[1];
sx q[1];
rz(-1.4453759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.611269) q[0];
sx q[0];
rz(-1.230254) q[0];
sx q[0];
rz(-2.5304619) q[0];
rz(-pi) q[1];
rz(-2.9640182) q[2];
sx q[2];
rz(-1.53731) q[2];
sx q[2];
rz(-2.356503) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0398324) q[1];
sx q[1];
rz(-1.2451485) q[1];
sx q[1];
rz(1.7937307) q[1];
rz(-1.2054005) q[3];
sx q[3];
rz(-0.8626079) q[3];
sx q[3];
rz(-2.303249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33862904) q[2];
sx q[2];
rz(-0.45946071) q[2];
sx q[2];
rz(-2.6098693) q[2];
rz(-1.7470597) q[3];
sx q[3];
rz(-1.8683542) q[3];
sx q[3];
rz(-1.1751706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8748473) q[0];
sx q[0];
rz(-1.5044455) q[0];
sx q[0];
rz(2.3496085) q[0];
rz(-2.2199471) q[1];
sx q[1];
rz(-1.6519203) q[1];
sx q[1];
rz(-2.0335061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1658679) q[0];
sx q[0];
rz(-0.97318968) q[0];
sx q[0];
rz(0.13522603) q[0];
rz(-pi) q[1];
rz(-3.0929977) q[2];
sx q[2];
rz(-0.94779769) q[2];
sx q[2];
rz(2.5041265) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82906065) q[1];
sx q[1];
rz(-0.8735114) q[1];
sx q[1];
rz(0.50364772) q[1];
rz(-1.4008303) q[3];
sx q[3];
rz(-0.32073354) q[3];
sx q[3];
rz(2.6331944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38829982) q[2];
sx q[2];
rz(-2.1123977) q[2];
sx q[2];
rz(1.8886867) q[2];
rz(-0.62475723) q[3];
sx q[3];
rz(-1.8936936) q[3];
sx q[3];
rz(-2.6784082) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6817634) q[0];
sx q[0];
rz(-2.9637931) q[0];
sx q[0];
rz(-0.27339992) q[0];
rz(-0.071648486) q[1];
sx q[1];
rz(-1.1337846) q[1];
sx q[1];
rz(1.515548) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26165379) q[0];
sx q[0];
rz(-2.0802705) q[0];
sx q[0];
rz(-1.1306056) q[0];
rz(1.927761) q[2];
sx q[2];
rz(-1.1276111) q[2];
sx q[2];
rz(0.24493327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7715986) q[1];
sx q[1];
rz(-2.150476) q[1];
sx q[1];
rz(1.7607776) q[1];
x q[2];
rz(1.203556) q[3];
sx q[3];
rz(-0.62588309) q[3];
sx q[3];
rz(-2.3123031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7800954) q[2];
sx q[2];
rz(-0.7146892) q[2];
sx q[2];
rz(-1.7717465) q[2];
rz(-1.0684377) q[3];
sx q[3];
rz(-1.3911824) q[3];
sx q[3];
rz(-2.1865602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84900981) q[0];
sx q[0];
rz(-0.42790258) q[0];
sx q[0];
rz(-3.0191315) q[0];
rz(0.017008688) q[1];
sx q[1];
rz(-1.5127134) q[1];
sx q[1];
rz(0.88517991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53986076) q[0];
sx q[0];
rz(-2.0575876) q[0];
sx q[0];
rz(2.9877325) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1954284) q[2];
sx q[2];
rz(-2.3705707) q[2];
sx q[2];
rz(-2.8147402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1840813) q[1];
sx q[1];
rz(-2.5169417) q[1];
sx q[1];
rz(-2.2923325) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2284352) q[3];
sx q[3];
rz(-0.71492787) q[3];
sx q[3];
rz(0.19107669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7077606) q[2];
sx q[2];
rz(-1.8261352) q[2];
sx q[2];
rz(3.0641277) q[2];
rz(2.7205983) q[3];
sx q[3];
rz(-1.9546031) q[3];
sx q[3];
rz(0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0635327) q[0];
sx q[0];
rz(-0.01883004) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(2.6497427) q[1];
sx q[1];
rz(-2.1147155) q[1];
sx q[1];
rz(1.9815725) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9710142) q[0];
sx q[0];
rz(-1.3828297) q[0];
sx q[0];
rz(2.2925799) q[0];
rz(2.5712058) q[2];
sx q[2];
rz(-1.6183934) q[2];
sx q[2];
rz(-0.80709761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8012538) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(-2.1992963) q[1];
rz(2.8468708) q[3];
sx q[3];
rz(-1.975276) q[3];
sx q[3];
rz(-1.3724788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.457095) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(2.1333466) q[2];
rz(1.6393939) q[3];
sx q[3];
rz(-1.1049756) q[3];
sx q[3];
rz(0.26386279) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098792583) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(-2.7368326) q[0];
rz(-0.81819355) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(-0.71664804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.19159) q[0];
sx q[0];
rz(-1.5806338) q[0];
sx q[0];
rz(-2.8623926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0990853) q[2];
sx q[2];
rz(-2.2620438) q[2];
sx q[2];
rz(-1.7886358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67988746) q[1];
sx q[1];
rz(-2.6033834) q[1];
sx q[1];
rz(2.6160282) q[1];
rz(-pi) q[2];
rz(-2.6022609) q[3];
sx q[3];
rz(-1.8908653) q[3];
sx q[3];
rz(2.0254997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8019668) q[2];
sx q[2];
rz(-2.5012987) q[2];
sx q[2];
rz(-0.65529811) q[2];
rz(2.7845434) q[3];
sx q[3];
rz(-1.7506295) q[3];
sx q[3];
rz(0.51957875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.201467) q[0];
sx q[0];
rz(-0.31620142) q[0];
sx q[0];
rz(1.0200208) q[0];
rz(-2.5992498) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(1.7592336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4750299) q[0];
sx q[0];
rz(-0.90797601) q[0];
sx q[0];
rz(2.3549838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21127659) q[2];
sx q[2];
rz(-1.4370434) q[2];
sx q[2];
rz(-2.0300031) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51782896) q[1];
sx q[1];
rz(-1.6490855) q[1];
sx q[1];
rz(2.0324367) q[1];
rz(3.0579733) q[3];
sx q[3];
rz(-1.5685387) q[3];
sx q[3];
rz(-1.797054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.072784) q[2];
sx q[2];
rz(-1.6135608) q[2];
sx q[2];
rz(-0.35045785) q[2];
rz(-0.46640486) q[3];
sx q[3];
rz(-0.91752183) q[3];
sx q[3];
rz(-0.94356999) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35895178) q[0];
sx q[0];
rz(-2.751001) q[0];
sx q[0];
rz(-2.1851831) q[0];
rz(1.6920754) q[1];
sx q[1];
rz(-0.78069514) q[1];
sx q[1];
rz(-0.67109674) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.886837) q[0];
sx q[0];
rz(-0.16567812) q[0];
sx q[0];
rz(1.0266233) q[0];
x q[1];
rz(-0.568957) q[2];
sx q[2];
rz(-1.7486186) q[2];
sx q[2];
rz(0.59019719) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8789323) q[1];
sx q[1];
rz(-1.9237681) q[1];
sx q[1];
rz(-3.073602) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70399474) q[3];
sx q[3];
rz(-1.3823095) q[3];
sx q[3];
rz(-2.2437167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87390071) q[2];
sx q[2];
rz(-0.42517069) q[2];
sx q[2];
rz(-2.0153866) q[2];
rz(2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(-3.0245074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0432334) q[0];
sx q[0];
rz(-2.5625304) q[0];
sx q[0];
rz(3.0845355) q[0];
rz(3.0077899) q[1];
sx q[1];
rz(-1.0452784) q[1];
sx q[1];
rz(0.71279508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063449115) q[0];
sx q[0];
rz(-2.8196555) q[0];
sx q[0];
rz(-0.19970317) q[0];
rz(-pi) q[1];
rz(2.7370896) q[2];
sx q[2];
rz(-1.9484011) q[2];
sx q[2];
rz(-1.9775569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9660898) q[1];
sx q[1];
rz(-1.1017297) q[1];
sx q[1];
rz(1.1203946) q[1];
rz(-pi) q[2];
rz(-3.0607575) q[3];
sx q[3];
rz(-2.2634215) q[3];
sx q[3];
rz(-2.5580542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5130676) q[2];
sx q[2];
rz(-1.8551989) q[2];
sx q[2];
rz(0.081550278) q[2];
rz(-1.9214123) q[3];
sx q[3];
rz(-0.27118513) q[3];
sx q[3];
rz(2.0512106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9995025) q[0];
sx q[0];
rz(-1.5220078) q[0];
sx q[0];
rz(-1.5600486) q[0];
rz(1.1355431) q[1];
sx q[1];
rz(-2.0027497) q[1];
sx q[1];
rz(-1.9948237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0479779) q[0];
sx q[0];
rz(-2.2824557) q[0];
sx q[0];
rz(-2.467903) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3421971) q[2];
sx q[2];
rz(-0.96916064) q[2];
sx q[2];
rz(0.55527273) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9830977) q[1];
sx q[1];
rz(-1.6545516) q[1];
sx q[1];
rz(0.2106481) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0542946) q[3];
sx q[3];
rz(-1.0074769) q[3];
sx q[3];
rz(2.6642852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5494988) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(0.36592323) q[2];
rz(-2.7538815) q[3];
sx q[3];
rz(-2.0230899) q[3];
sx q[3];
rz(1.6754735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6112919) q[0];
sx q[0];
rz(-0.74321754) q[0];
sx q[0];
rz(2.8331953) q[0];
rz(-0.74956924) q[1];
sx q[1];
rz(-1.1309962) q[1];
sx q[1];
rz(1.8684734) q[1];
rz(0.91084008) q[2];
sx q[2];
rz(-1.2963933) q[2];
sx q[2];
rz(-2.8252841) q[2];
rz(1.9182792) q[3];
sx q[3];
rz(-0.96746222) q[3];
sx q[3];
rz(-0.37012561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
