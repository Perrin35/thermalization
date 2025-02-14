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
rz(-0.24124423) q[0];
sx q[0];
rz(-0.21114199) q[0];
sx q[0];
rz(0.40786064) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(3.0768375) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252547) q[0];
sx q[0];
rz(-3.0260575) q[0];
sx q[0];
rz(-1.6166685) q[0];
x q[1];
rz(-0.98528905) q[2];
sx q[2];
rz(-1.3006214) q[2];
sx q[2];
rz(-0.5952685) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25542361) q[1];
sx q[1];
rz(-1.5597289) q[1];
sx q[1];
rz(-0.14126417) q[1];
rz(-2.519386) q[3];
sx q[3];
rz(-1.6554621) q[3];
sx q[3];
rz(-1.0877999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9593418) q[2];
sx q[2];
rz(-1.094341) q[2];
sx q[2];
rz(2.0056637) q[2];
rz(2.411339) q[3];
sx q[3];
rz(-1.3948995) q[3];
sx q[3];
rz(1.6212538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0486883) q[0];
sx q[0];
rz(-0.95656675) q[0];
sx q[0];
rz(-3.0862869) q[0];
rz(-0.37120184) q[1];
sx q[1];
rz(-2.7822045) q[1];
sx q[1];
rz(2.7296383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13875554) q[0];
sx q[0];
rz(-1.5056207) q[0];
sx q[0];
rz(0.52703339) q[0];
rz(1.2367593) q[2];
sx q[2];
rz(-1.8530117) q[2];
sx q[2];
rz(1.7441234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.72295475) q[1];
sx q[1];
rz(-1.4297863) q[1];
sx q[1];
rz(-2.4055071) q[1];
rz(-2.7633414) q[3];
sx q[3];
rz(-2.1839199) q[3];
sx q[3];
rz(-0.65968266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61341316) q[2];
sx q[2];
rz(-1.6772905) q[2];
sx q[2];
rz(0.8826274) q[2];
rz(-1.708185) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(-0.17361704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16647896) q[0];
sx q[0];
rz(-3.1298895) q[0];
sx q[0];
rz(0.56667462) q[0];
rz(-2.7239679) q[1];
sx q[1];
rz(-1.5387225) q[1];
sx q[1];
rz(1.3272939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8493932) q[0];
sx q[0];
rz(-2.4367094) q[0];
sx q[0];
rz(-1.7017176) q[0];
rz(2.9887373) q[2];
sx q[2];
rz(-2.1899583) q[2];
sx q[2];
rz(-2.4094606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.23073634) q[1];
sx q[1];
rz(-1.5720782) q[1];
sx q[1];
rz(1.8788473) q[1];
rz(2.684115) q[3];
sx q[3];
rz(-0.8972392) q[3];
sx q[3];
rz(2.9411773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7495482) q[2];
sx q[2];
rz(-2.0508524) q[2];
sx q[2];
rz(-2.2331179) q[2];
rz(2.3467017) q[3];
sx q[3];
rz(-1.1029693) q[3];
sx q[3];
rz(-0.046549646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6598776) q[0];
sx q[0];
rz(-0.91232038) q[0];
sx q[0];
rz(0.6657486) q[0];
rz(2.2944229) q[1];
sx q[1];
rz(-0.24996346) q[1];
sx q[1];
rz(-0.4096823) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2862662) q[0];
sx q[0];
rz(-0.5200035) q[0];
sx q[0];
rz(-2.8306243) q[0];
x q[1];
rz(1.9471859) q[2];
sx q[2];
rz(-2.3431398) q[2];
sx q[2];
rz(1.2898766) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53654251) q[1];
sx q[1];
rz(-2.4808421) q[1];
sx q[1];
rz(-1.9416757) q[1];
rz(1.5104896) q[3];
sx q[3];
rz(-2.6968287) q[3];
sx q[3];
rz(0.30032762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72982558) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(-1.3789619) q[2];
rz(-1.1369368) q[3];
sx q[3];
rz(-1.3552908) q[3];
sx q[3];
rz(2.304145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9627422) q[0];
sx q[0];
rz(-1.0609635) q[0];
sx q[0];
rz(-0.41819292) q[0];
rz(1.4232945) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(-1.4994015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9327576) q[0];
sx q[0];
rz(-1.3262188) q[0];
sx q[0];
rz(0.45758684) q[0];
rz(-pi) q[1];
rz(-0.73694456) q[2];
sx q[2];
rz(-0.43790753) q[2];
sx q[2];
rz(-2.2543668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.114257) q[1];
sx q[1];
rz(-1.3523977) q[1];
sx q[1];
rz(-0.059192358) q[1];
x q[2];
rz(-3.0762227) q[3];
sx q[3];
rz(-1.7002859) q[3];
sx q[3];
rz(-2.0355899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2749918) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(-0.17821136) q[2];
rz(2.7675659) q[3];
sx q[3];
rz(-2.1686797) q[3];
sx q[3];
rz(1.1427574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571534) q[0];
sx q[0];
rz(-1.8171808) q[0];
sx q[0];
rz(-1.4763747) q[0];
rz(-2.5560675) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(-2.4077328) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7134345) q[0];
sx q[0];
rz(-2.1913986) q[0];
sx q[0];
rz(-2.1290859) q[0];
x q[1];
rz(-0.30690212) q[2];
sx q[2];
rz(-1.2595121) q[2];
sx q[2];
rz(2.5800623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3702995) q[1];
sx q[1];
rz(-1.7669189) q[1];
sx q[1];
rz(2.0972154) q[1];
rz(-pi) q[2];
rz(2.0229983) q[3];
sx q[3];
rz(-1.7590932) q[3];
sx q[3];
rz(1.888243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0840941) q[2];
sx q[2];
rz(-1.3749264) q[2];
sx q[2];
rz(1.6424087) q[2];
rz(1.621014) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(-2.9852941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3868189) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(-1.3936438) q[0];
rz(-2.5539894) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(0.29744068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818243) q[0];
sx q[0];
rz(-2.3399669) q[0];
sx q[0];
rz(-1.0645435) q[0];
rz(-pi) q[1];
rz(-1.5027375) q[2];
sx q[2];
rz(-0.27114332) q[2];
sx q[2];
rz(-0.23672297) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.12402) q[1];
sx q[1];
rz(-2.4618853) q[1];
sx q[1];
rz(1.555383) q[1];
rz(-pi) q[2];
rz(1.8930412) q[3];
sx q[3];
rz(-0.6686223) q[3];
sx q[3];
rz(3.0361255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3191159) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(-2.759867) q[2];
rz(2.5413359) q[3];
sx q[3];
rz(-0.67265066) q[3];
sx q[3];
rz(0.265358) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6622019) q[0];
sx q[0];
rz(-1.6265765) q[0];
sx q[0];
rz(-1.8040682) q[0];
rz(-3.1321101) q[1];
sx q[1];
rz(-2.1818706) q[1];
sx q[1];
rz(1.7705852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66390304) q[0];
sx q[0];
rz(-0.36720095) q[0];
sx q[0];
rz(2.7504671) q[0];
rz(-0.96641936) q[2];
sx q[2];
rz(-0.76497173) q[2];
sx q[2];
rz(1.6993831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7035455) q[1];
sx q[1];
rz(-2.3444972) q[1];
sx q[1];
rz(0.17532562) q[1];
rz(-pi) q[2];
rz(1.480732) q[3];
sx q[3];
rz(-1.1582631) q[3];
sx q[3];
rz(1.5145258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3235772) q[2];
sx q[2];
rz(-2.24021) q[2];
sx q[2];
rz(0.14249194) q[2];
rz(-2.2583708) q[3];
sx q[3];
rz(-1.764069) q[3];
sx q[3];
rz(-1.4628598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.091081) q[0];
sx q[0];
rz(-0.091878042) q[0];
sx q[0];
rz(-1.2485414) q[0];
rz(2.7342791) q[1];
sx q[1];
rz(-1.3469603) q[1];
sx q[1];
rz(1.1936845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68797217) q[0];
sx q[0];
rz(-1.2415452) q[0];
sx q[0];
rz(-2.8696174) q[0];
x q[1];
rz(0.87393729) q[2];
sx q[2];
rz(-0.94183445) q[2];
sx q[2];
rz(-3.0203117) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9858169) q[1];
sx q[1];
rz(-1.078842) q[1];
sx q[1];
rz(-1.953275) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3109929) q[3];
sx q[3];
rz(-1.570829) q[3];
sx q[3];
rz(-2.3268229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.074389) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(1.1129145) q[2];
rz(-3.0985966) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-0.23785166) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766895) q[0];
sx q[0];
rz(-1.4142798) q[0];
sx q[0];
rz(-0.195737) q[0];
rz(1.9211357) q[1];
sx q[1];
rz(-2.7595322) q[1];
sx q[1];
rz(-2.1641796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5479056) q[0];
sx q[0];
rz(-1.4713877) q[0];
sx q[0];
rz(0.25940827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1117729) q[2];
sx q[2];
rz(-2.8089097) q[2];
sx q[2];
rz(-2.7667342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43454506) q[1];
sx q[1];
rz(-2.0710398) q[1];
sx q[1];
rz(-2.2342199) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39267003) q[3];
sx q[3];
rz(-1.8373946) q[3];
sx q[3];
rz(-2.6368898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0105373) q[2];
sx q[2];
rz(-0.26630339) q[2];
sx q[2];
rz(-1.4492501) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-1.0613469) q[3];
sx q[3];
rz(-1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5129678) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.3799008) q[1];
sx q[1];
rz(-1.6086144) q[1];
sx q[1];
rz(-0.10133941) q[1];
rz(-2.9607282) q[2];
sx q[2];
rz(-1.7500063) q[2];
sx q[2];
rz(1.0371006) q[2];
rz(-0.32573777) q[3];
sx q[3];
rz(-0.67254638) q[3];
sx q[3];
rz(-1.2253958) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
