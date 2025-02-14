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
rz(2.7316982) q[0];
sx q[0];
rz(-1.7562261) q[0];
sx q[0];
rz(2.07055) q[0];
rz(-1.3261803) q[1];
sx q[1];
rz(-2.2161127) q[1];
sx q[1];
rz(-2.7993536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5919094) q[0];
sx q[0];
rz(-0.91381493) q[0];
sx q[0];
rz(0.4763989) q[0];
rz(2.3314807) q[2];
sx q[2];
rz(-1.9829921) q[2];
sx q[2];
rz(-0.10998943) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46821478) q[1];
sx q[1];
rz(-1.4644551) q[1];
sx q[1];
rz(-2.0061226) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5299523) q[3];
sx q[3];
rz(-0.39642912) q[3];
sx q[3];
rz(1.2547697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60275117) q[2];
sx q[2];
rz(-2.4138236) q[2];
sx q[2];
rz(-1.4347428) q[2];
rz(2.1753066) q[3];
sx q[3];
rz(-1.9478226) q[3];
sx q[3];
rz(-0.57893354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.12427881) q[0];
sx q[0];
rz(-2.3983428) q[0];
sx q[0];
rz(-3.0857575) q[0];
rz(2.3827379) q[1];
sx q[1];
rz(-1.1412303) q[1];
sx q[1];
rz(0.32078823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082328) q[0];
sx q[0];
rz(-0.77455168) q[0];
sx q[0];
rz(-2.9288128) q[0];
rz(-2.3622038) q[2];
sx q[2];
rz(-1.3192799) q[2];
sx q[2];
rz(0.0096461065) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4652511) q[1];
sx q[1];
rz(-0.23561978) q[1];
sx q[1];
rz(1.444677) q[1];
x q[2];
rz(2.0386215) q[3];
sx q[3];
rz(-2.247034) q[3];
sx q[3];
rz(-2.6060864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8212829) q[2];
sx q[2];
rz(-1.8822957) q[2];
sx q[2];
rz(-0.28309509) q[2];
rz(0.29575944) q[3];
sx q[3];
rz(-1.4247954) q[3];
sx q[3];
rz(-1.504771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54842424) q[0];
sx q[0];
rz(-1.230509) q[0];
sx q[0];
rz(3.044627) q[0];
rz(-1.7123669) q[1];
sx q[1];
rz(-1.2173419) q[1];
sx q[1];
rz(2.7860577) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9331207) q[0];
sx q[0];
rz(-1.2829395) q[0];
sx q[0];
rz(-2.5639749) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41496241) q[2];
sx q[2];
rz(-1.6646766) q[2];
sx q[2];
rz(1.0549269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.695836) q[1];
sx q[1];
rz(-1.9805659) q[1];
sx q[1];
rz(-0.10182937) q[1];
rz(-pi) q[2];
rz(2.8899651) q[3];
sx q[3];
rz(-1.5306002) q[3];
sx q[3];
rz(-0.1110368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43704438) q[2];
sx q[2];
rz(-1.9883678) q[2];
sx q[2];
rz(-2.8538749) q[2];
rz(-0.25117609) q[3];
sx q[3];
rz(-2.0468678) q[3];
sx q[3];
rz(-1.8672966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1890202) q[0];
sx q[0];
rz(-2.548521) q[0];
sx q[0];
rz(0.75743842) q[0];
rz(2.3144552) q[1];
sx q[1];
rz(-0.65991455) q[1];
sx q[1];
rz(-1.3077024) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33117244) q[0];
sx q[0];
rz(-3.0068827) q[0];
sx q[0];
rz(-1.170838) q[0];
rz(2.1266379) q[2];
sx q[2];
rz(-0.53947811) q[2];
sx q[2];
rz(-1.2459823) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2590307) q[1];
sx q[1];
rz(-0.81444959) q[1];
sx q[1];
rz(-2.7689054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3263474) q[3];
sx q[3];
rz(-1.5119879) q[3];
sx q[3];
rz(1.5016426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7149675) q[2];
sx q[2];
rz(-1.7968618) q[2];
sx q[2];
rz(-1.7460543) q[2];
rz(1.7914145) q[3];
sx q[3];
rz(-1.502864) q[3];
sx q[3];
rz(-3.1135528) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79513079) q[0];
sx q[0];
rz(-0.53855723) q[0];
sx q[0];
rz(1.7002456) q[0];
rz(-2.2137008) q[1];
sx q[1];
rz(-0.83506942) q[1];
sx q[1];
rz(2.3477614) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071811) q[0];
sx q[0];
rz(-0.4259122) q[0];
sx q[0];
rz(-1.3176729) q[0];
rz(-1.687811) q[2];
sx q[2];
rz(-1.4037366) q[2];
sx q[2];
rz(-2.2599205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9873569) q[1];
sx q[1];
rz(-1.8786228) q[1];
sx q[1];
rz(-1.0197784) q[1];
x q[2];
rz(-0.63266261) q[3];
sx q[3];
rz(-2.5404937) q[3];
sx q[3];
rz(-0.68279235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64357197) q[2];
sx q[2];
rz(-1.6732432) q[2];
sx q[2];
rz(0.35430655) q[2];
rz(-2.0453359) q[3];
sx q[3];
rz(-2.8772964) q[3];
sx q[3];
rz(-1.8335584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7588014) q[0];
sx q[0];
rz(-2.379874) q[0];
sx q[0];
rz(0.85269165) q[0];
rz(2.8755152) q[1];
sx q[1];
rz(-2.1177025) q[1];
sx q[1];
rz(1.783225) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.81937) q[0];
sx q[0];
rz(-2.529105) q[0];
sx q[0];
rz(1.9507381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9044962) q[2];
sx q[2];
rz(-0.70095567) q[2];
sx q[2];
rz(-2.440941) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0222802) q[1];
sx q[1];
rz(-1.90442) q[1];
sx q[1];
rz(-2.5319478) q[1];
rz(-pi) q[2];
rz(0.27652506) q[3];
sx q[3];
rz(-0.92036906) q[3];
sx q[3];
rz(2.1610633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69409662) q[2];
sx q[2];
rz(-1.2022377) q[2];
sx q[2];
rz(2.6608432) q[2];
rz(-2.2065744) q[3];
sx q[3];
rz(-2.1604249) q[3];
sx q[3];
rz(2.7267314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9347436) q[0];
sx q[0];
rz(-0.54731363) q[0];
sx q[0];
rz(-0.68688399) q[0];
rz(0.96784776) q[1];
sx q[1];
rz(-1.5279488) q[1];
sx q[1];
rz(0.17766775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8657836) q[0];
sx q[0];
rz(-1.1622419) q[0];
sx q[0];
rz(-1.1955117) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91089852) q[2];
sx q[2];
rz(-1.0086806) q[2];
sx q[2];
rz(-2.0695994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2100169) q[1];
sx q[1];
rz(-2.615395) q[1];
sx q[1];
rz(0.99858649) q[1];
rz(-pi) q[2];
rz(1.6611886) q[3];
sx q[3];
rz(-1.6393344) q[3];
sx q[3];
rz(-1.3892421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.054242) q[2];
sx q[2];
rz(-2.6998616) q[2];
sx q[2];
rz(2.1596215) q[2];
rz(0.25518498) q[3];
sx q[3];
rz(-1.988215) q[3];
sx q[3];
rz(2.065778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0690689) q[0];
sx q[0];
rz(-2.531932) q[0];
sx q[0];
rz(0.12741086) q[0];
rz(1.472507) q[1];
sx q[1];
rz(-1.1839048) q[1];
sx q[1];
rz(-1.4471819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54501891) q[0];
sx q[0];
rz(-1.8880287) q[0];
sx q[0];
rz(-1.7659643) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54355343) q[2];
sx q[2];
rz(-1.2464036) q[2];
sx q[2];
rz(-2.3860562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.060272) q[1];
sx q[1];
rz(-1.6439207) q[1];
sx q[1];
rz(3.0036219) q[1];
rz(-pi) q[2];
x q[2];
rz(2.892268) q[3];
sx q[3];
rz(-2.1405947) q[3];
sx q[3];
rz(2.9227481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.24826) q[2];
sx q[2];
rz(-2.1620763) q[2];
sx q[2];
rz(2.7323006) q[2];
rz(0.69581699) q[3];
sx q[3];
rz(-2.4455363) q[3];
sx q[3];
rz(2.5223562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61243397) q[0];
sx q[0];
rz(-2.9127064) q[0];
sx q[0];
rz(2.9869475) q[0];
rz(2.0507428) q[1];
sx q[1];
rz(-2.2583074) q[1];
sx q[1];
rz(-1.326391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94448254) q[0];
sx q[0];
rz(-1.3061227) q[0];
sx q[0];
rz(-0.84828429) q[0];
rz(-pi) q[1];
rz(-3.0934264) q[2];
sx q[2];
rz(-1.0720952) q[2];
sx q[2];
rz(1.1737385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.058097046) q[1];
sx q[1];
rz(-2.3196583) q[1];
sx q[1];
rz(0.067598299) q[1];
x q[2];
rz(-0.38549586) q[3];
sx q[3];
rz(-0.49043819) q[3];
sx q[3];
rz(2.3602594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0632625) q[2];
sx q[2];
rz(-0.41389725) q[2];
sx q[2];
rz(-2.4388893) q[2];
rz(-3.0472158) q[3];
sx q[3];
rz(-2.1636212) q[3];
sx q[3];
rz(0.92488658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6150045) q[0];
sx q[0];
rz(-2.1115392) q[0];
sx q[0];
rz(-2.7870542) q[0];
rz(-1.7189369) q[1];
sx q[1];
rz(-1.4965897) q[1];
sx q[1];
rz(2.6197701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3238597) q[0];
sx q[0];
rz(-2.7561152) q[0];
sx q[0];
rz(-1.9800703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9873677) q[2];
sx q[2];
rz(-2.8099001) q[2];
sx q[2];
rz(2.7046534) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21510151) q[1];
sx q[1];
rz(-2.1128078) q[1];
sx q[1];
rz(-0.78185977) q[1];
rz(-2.4168968) q[3];
sx q[3];
rz(-1.0130289) q[3];
sx q[3];
rz(-0.23797056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78174138) q[2];
sx q[2];
rz(-1.968911) q[2];
sx q[2];
rz(-2.4594128) q[2];
rz(1.773268) q[3];
sx q[3];
rz(-1.5409639) q[3];
sx q[3];
rz(-2.1413596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.901004) q[0];
sx q[0];
rz(-0.95745845) q[0];
sx q[0];
rz(-0.29314713) q[0];
rz(-2.3121569) q[1];
sx q[1];
rz(-1.4301626) q[1];
sx q[1];
rz(1.593874) q[1];
rz(-1.310036) q[2];
sx q[2];
rz(-1.3481067) q[2];
sx q[2];
rz(1.3630661) q[2];
rz(-2.8021723) q[3];
sx q[3];
rz(-0.11133736) q[3];
sx q[3];
rz(2.8471024) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
