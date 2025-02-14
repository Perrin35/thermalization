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
rz(3.1495659) q[0];
sx q[0];
rz(6.719448) q[0];
sx q[0];
rz(8.414581) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(7.0893256) q[1];
sx q[1];
rz(10.448699) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7069081) q[0];
sx q[0];
rz(-0.48792612) q[0];
sx q[0];
rz(1.5561681) q[0];
rz(-pi) q[1];
rz(-1.3390117) q[2];
sx q[2];
rz(-1.9531774) q[2];
sx q[2];
rz(0.2222375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6758821) q[1];
sx q[1];
rz(-2.9799068) q[1];
sx q[1];
rz(-1.449701) q[1];
x q[2];
rz(0.2497345) q[3];
sx q[3];
rz(-1.6563376) q[3];
sx q[3];
rz(-1.6289323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1930834) q[2];
sx q[2];
rz(-0.32863363) q[2];
sx q[2];
rz(0.23552093) q[2];
rz(-1.0282907) q[3];
sx q[3];
rz(-1.8833269) q[3];
sx q[3];
rz(-1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.704772) q[0];
sx q[0];
rz(-1.3587767) q[0];
sx q[0];
rz(0.92192465) q[0];
rz(-3.0251265) q[1];
sx q[1];
rz(-1.7200229) q[1];
sx q[1];
rz(2.2540653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5619316) q[0];
sx q[0];
rz(-1.4401739) q[0];
sx q[0];
rz(-1.7586673) q[0];
x q[1];
rz(2.9824184) q[2];
sx q[2];
rz(-2.7870057) q[2];
sx q[2];
rz(-1.872242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.973399) q[1];
sx q[1];
rz(-1.0032986) q[1];
sx q[1];
rz(-1.7048111) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55887115) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(2.9026046) q[2];
rz(2.41921) q[3];
sx q[3];
rz(-2.3014849) q[3];
sx q[3];
rz(1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8257985) q[0];
sx q[0];
rz(-0.9676942) q[0];
sx q[0];
rz(0.78829366) q[0];
rz(1.7514508) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(-1.4424666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628796) q[0];
sx q[0];
rz(-1.4294693) q[0];
sx q[0];
rz(3.0706329) q[0];
x q[1];
rz(2.2295775) q[2];
sx q[2];
rz(-0.46123966) q[2];
sx q[2];
rz(1.3818503) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6953814) q[1];
sx q[1];
rz(-1.808481) q[1];
sx q[1];
rz(-2.2878245) q[1];
rz(-1.0679108) q[3];
sx q[3];
rz(-1.4137795) q[3];
sx q[3];
rz(2.5105421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4060789) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(-0.25685143) q[2];
rz(0.86366051) q[3];
sx q[3];
rz(-1.0830027) q[3];
sx q[3];
rz(-2.5679722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470806) q[0];
sx q[0];
rz(-3.1086476) q[0];
sx q[0];
rz(-1.5091913) q[0];
rz(-0.31653658) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(-1.0155771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7596801) q[0];
sx q[0];
rz(-2.0354871) q[0];
sx q[0];
rz(-2.0813991) q[0];
rz(-2.8695753) q[2];
sx q[2];
rz(-0.81230703) q[2];
sx q[2];
rz(2.7886632) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9719661) q[1];
sx q[1];
rz(-1.7952654) q[1];
sx q[1];
rz(1.1826993) q[1];
rz(-pi) q[2];
rz(-0.16839858) q[3];
sx q[3];
rz(-2.3108498) q[3];
sx q[3];
rz(0.5357045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(-1.451937) q[2];
rz(1.9076094) q[3];
sx q[3];
rz(-1.0943639) q[3];
sx q[3];
rz(2.2598677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.7538309) q[0];
sx q[0];
rz(-1.8330638) q[0];
sx q[0];
rz(1.0863139) q[0];
rz(1.7408675) q[1];
sx q[1];
rz(-0.81175214) q[1];
sx q[1];
rz(-0.0033671826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2732179) q[0];
sx q[0];
rz(-1.5649319) q[0];
sx q[0];
rz(-0.59451367) q[0];
rz(-2.1589627) q[2];
sx q[2];
rz(-1.2890883) q[2];
sx q[2];
rz(-2.1717193) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7351384) q[1];
sx q[1];
rz(-0.5105831) q[1];
sx q[1];
rz(3.0264454) q[1];
rz(-0.037854511) q[3];
sx q[3];
rz(-2.541171) q[3];
sx q[3];
rz(0.63943938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9747666) q[2];
sx q[2];
rz(-0.70009118) q[2];
sx q[2];
rz(0.33494803) q[2];
rz(0.29609984) q[3];
sx q[3];
rz(-2.1985603) q[3];
sx q[3];
rz(3.0193442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.629338) q[0];
sx q[0];
rz(-0.4991931) q[0];
sx q[0];
rz(0.41803023) q[0];
rz(0.66608518) q[1];
sx q[1];
rz(-1.8981551) q[1];
sx q[1];
rz(-0.24806771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46588308) q[0];
sx q[0];
rz(-2.0572691) q[0];
sx q[0];
rz(-0.37635641) q[0];
rz(-pi) q[1];
rz(-1.4573967) q[2];
sx q[2];
rz(-0.45219401) q[2];
sx q[2];
rz(-0.8314641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0478749) q[1];
sx q[1];
rz(-1.1823726) q[1];
sx q[1];
rz(1.018143) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2002688) q[3];
sx q[3];
rz(-1.7793312) q[3];
sx q[3];
rz(2.1311614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.69998133) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(-1.1629533) q[2];
rz(0.59398389) q[3];
sx q[3];
rz(-1.6006399) q[3];
sx q[3];
rz(2.1514413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.537259) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(-0.064924031) q[0];
rz(-0.36554947) q[1];
sx q[1];
rz(-1.431501) q[1];
sx q[1];
rz(0.72986594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54016544) q[0];
sx q[0];
rz(-1.6807204) q[0];
sx q[0];
rz(-0.92782995) q[0];
rz(-pi) q[1];
rz(2.9403749) q[2];
sx q[2];
rz(-2.3683753) q[2];
sx q[2];
rz(1.4035483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0893475) q[1];
sx q[1];
rz(-0.56530732) q[1];
sx q[1];
rz(-1.7491399) q[1];
x q[2];
rz(-1.1564674) q[3];
sx q[3];
rz(-0.92770139) q[3];
sx q[3];
rz(-0.82655686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87021747) q[2];
sx q[2];
rz(-1.7729746) q[2];
sx q[2];
rz(-2.0591056) q[2];
rz(1.4879976) q[3];
sx q[3];
rz(-0.36181417) q[3];
sx q[3];
rz(-2.7910119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38700405) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(0.3311232) q[0];
rz(1.4777199) q[1];
sx q[1];
rz(-2.176087) q[1];
sx q[1];
rz(-2.8211735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77003682) q[0];
sx q[0];
rz(-1.5450243) q[0];
sx q[0];
rz(0.098866247) q[0];
rz(0.75866239) q[2];
sx q[2];
rz(-2.9305923) q[2];
sx q[2];
rz(1.3727486) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81605212) q[1];
sx q[1];
rz(-2.7150702) q[1];
sx q[1];
rz(-0.037530516) q[1];
rz(1.5005169) q[3];
sx q[3];
rz(-0.68551862) q[3];
sx q[3];
rz(-2.7932515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7122345) q[2];
sx q[2];
rz(-1.4171436) q[2];
sx q[2];
rz(-2.7799535) q[2];
rz(-2.2278191) q[3];
sx q[3];
rz(-2.8445966) q[3];
sx q[3];
rz(2.698212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2509895) q[0];
sx q[0];
rz(-2.3539703) q[0];
sx q[0];
rz(0.11601624) q[0];
rz(1.3277671) q[1];
sx q[1];
rz(-1.1632183) q[1];
sx q[1];
rz(-1.2001002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2105149) q[0];
sx q[0];
rz(-1.7621982) q[0];
sx q[0];
rz(-1.5516419) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88516) q[2];
sx q[2];
rz(-1.6994393) q[2];
sx q[2];
rz(0.41510116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1397651) q[1];
sx q[1];
rz(-1.5561625) q[1];
sx q[1];
rz(0.56165265) q[1];
x q[2];
rz(-0.16606826) q[3];
sx q[3];
rz(-1.8072053) q[3];
sx q[3];
rz(2.0061559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8414595) q[2];
sx q[2];
rz(-1.4347142) q[2];
sx q[2];
rz(0.33162281) q[2];
rz(-2.8578109) q[3];
sx q[3];
rz(-1.1016568) q[3];
sx q[3];
rz(-0.42154977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83061853) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(0.16657883) q[0];
rz(2.2919948) q[1];
sx q[1];
rz(-0.46509898) q[1];
sx q[1];
rz(0.073624484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62944996) q[0];
sx q[0];
rz(-2.3733022) q[0];
sx q[0];
rz(0.37017249) q[0];
rz(-pi) q[1];
rz(0.4771941) q[2];
sx q[2];
rz(-1.0721803) q[2];
sx q[2];
rz(2.3248364) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4576118) q[1];
sx q[1];
rz(-1.906412) q[1];
sx q[1];
rz(-0.5313576) q[1];
rz(-pi) q[2];
rz(-1.0898548) q[3];
sx q[3];
rz(-2.1409799) q[3];
sx q[3];
rz(2.7358304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2354043) q[2];
sx q[2];
rz(-2.3222458) q[2];
sx q[2];
rz(-2.4833615) q[2];
rz(-1.341358) q[3];
sx q[3];
rz(-2.1737289) q[3];
sx q[3];
rz(-0.85097504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8258719) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(0.52147621) q[1];
sx q[1];
rz(-1.5745402) q[1];
sx q[1];
rz(1.5706617) q[1];
rz(-2.3519631) q[2];
sx q[2];
rz(-1.9676002) q[2];
sx q[2];
rz(-1.0330539) q[2];
rz(2.519532) q[3];
sx q[3];
rz(-1.4348772) q[3];
sx q[3];
rz(0.025195599) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
