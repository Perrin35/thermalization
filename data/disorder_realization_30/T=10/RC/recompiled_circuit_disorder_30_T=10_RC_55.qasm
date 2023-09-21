OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49139872) q[0];
sx q[0];
rz(-0.2645275) q[0];
sx q[0];
rz(-0.39443031) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(2.8013464) q[1];
sx q[1];
rz(10.624788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5651363) q[0];
sx q[0];
rz(1.5229043) q[0];
x q[1];
rz(2.8067402) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(1.2644757) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29772273) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-0.86617275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11043926) q[3];
sx q[3];
rz(-2.7613598) q[3];
sx q[3];
rz(1.5766174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(0.34399024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-0.66939676) q[1];
sx q[1];
rz(1.7864236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38793135) q[0];
sx q[0];
rz(-2.5270215) q[0];
sx q[0];
rz(-1.4580926) q[0];
x q[1];
rz(-1.4488892) q[2];
sx q[2];
rz(-1.4093471) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.018108798) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(2.0204087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7543206) q[3];
sx q[3];
rz(-1.6471383) q[3];
sx q[3];
rz(2.9294088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(2.6039092) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(2.5007201) q[0];
rz(-0.74869853) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(2.0764988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5168034) q[0];
sx q[0];
rz(-1.5644801) q[0];
sx q[0];
rz(1.9301231) q[0];
rz(-pi) q[1];
rz(-2.2674019) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(0.092560571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.47485891) q[1];
sx q[1];
rz(-1.6068646) q[1];
sx q[1];
rz(-0.95400793) q[1];
rz(-0.25032708) q[3];
sx q[3];
rz(-0.86391376) q[3];
sx q[3];
rz(1.3644497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.179402) q[2];
sx q[2];
rz(-0.71737814) q[2];
rz(-0.68850368) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(2.8919019) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(-3.1304741) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5866833) q[0];
sx q[0];
rz(-2.0218098) q[0];
sx q[0];
rz(-2.6888072) q[0];
rz(-pi) q[1];
x q[1];
rz(2.501802) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(-1.7196136) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3938155) q[1];
sx q[1];
rz(-1.2443466) q[1];
sx q[1];
rz(2.087489) q[1];
x q[2];
rz(-0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(-1.1312248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6461688) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(-2.8584976) q[2];
rz(-2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-2.6470673) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.8146851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4715695) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(0.75385401) q[0];
rz(-pi) q[1];
rz(2.3438498) q[2];
sx q[2];
rz(-2.0125055) q[2];
sx q[2];
rz(-1.1898578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27969589) q[1];
sx q[1];
rz(-1.3394636) q[1];
sx q[1];
rz(0.33613236) q[1];
rz(0.29019659) q[3];
sx q[3];
rz(-1.962933) q[3];
sx q[3];
rz(0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.2456606) q[2];
rz(-1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(0.6814878) q[0];
rz(-2.9340414) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(2.1372165) q[0];
rz(-pi) q[1];
rz(2.240928) q[2];
sx q[2];
rz(-1.2256983) q[2];
sx q[2];
rz(-0.64054856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7129732) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(0.58847217) q[1];
x q[2];
rz(-2.9168105) q[3];
sx q[3];
rz(-2.4596679) q[3];
sx q[3];
rz(-0.080760591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(2.7872655) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(-2.7697146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-2.4672467) q[0];
rz(-2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(0.56232125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494941) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(2.2390319) q[0];
rz(-pi) q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-2.7966768) q[2];
sx q[2];
rz(3.0722741) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79531419) q[1];
sx q[1];
rz(-0.242713) q[1];
sx q[1];
rz(0.64530356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.167422) q[3];
sx q[3];
rz(-1.4985634) q[3];
sx q[3];
rz(0.40856397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0397296) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(-0.34004655) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(-3.0026644) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.6202392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28946653) q[0];
sx q[0];
rz(-0.20848256) q[0];
sx q[0];
rz(2.8058488) q[0];
x q[1];
rz(2.0422158) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(-0.21381703) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.382114) q[1];
sx q[1];
rz(-1.5112875) q[1];
sx q[1];
rz(0.77854034) q[1];
rz(-pi) q[2];
rz(-2.2500854) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(-2.8241983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-0.29433027) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.3785988) q[3];
sx q[3];
rz(0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(0.35933581) q[0];
rz(2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(-2.8709581) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198869) q[0];
sx q[0];
rz(-1.5396376) q[0];
sx q[0];
rz(-3.1057538) q[0];
rz(-pi) q[1];
rz(2.3324899) q[2];
sx q[2];
rz(-1.9168233) q[2];
sx q[2];
rz(-1.5835294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.86079019) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(1.6153107) q[1];
rz(0.6587894) q[3];
sx q[3];
rz(-2.3156392) q[3];
sx q[3];
rz(-2.7487019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3140807) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(-0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(2.5922095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.6332508) q[0];
sx q[0];
rz(-2.4023138) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(2.646692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65942818) q[0];
sx q[0];
rz(-2.6539301) q[0];
sx q[0];
rz(1.6517261) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8137389) q[2];
sx q[2];
rz(-1.5640537) q[2];
sx q[2];
rz(1.1525796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1315688) q[1];
sx q[1];
rz(-1.459517) q[1];
sx q[1];
rz(-2.2113423) q[1];
x q[2];
rz(-1.4570191) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(-2.8137394) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(-1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(2.1140425) q[2];
sx q[2];
rz(-1.7474568) q[2];
sx q[2];
rz(-2.2251868) q[2];
rz(-0.26631793) q[3];
sx q[3];
rz(-1.3238293) q[3];
sx q[3];
rz(-2.7662591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
