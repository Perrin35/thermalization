OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(2.7881665) q[0];
sx q[0];
rz(8.3600144) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(-0.53607166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85815) q[0];
sx q[0];
rz(-1.802823) q[0];
sx q[0];
rz(2.0767077) q[0];
rz(-pi) q[1];
rz(-0.86059086) q[2];
sx q[2];
rz(-1.6221559) q[2];
sx q[2];
rz(-2.0757338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88692611) q[1];
sx q[1];
rz(-2.1517422) q[1];
sx q[1];
rz(1.3455774) q[1];
rz(-1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(-2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(2.8895767) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(2.305796) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(2.4332411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3591374) q[0];
sx q[0];
rz(-1.29175) q[0];
sx q[0];
rz(0.38048394) q[0];
rz(-2.9950036) q[2];
sx q[2];
rz(-0.97448889) q[2];
sx q[2];
rz(1.6080315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2409776) q[1];
sx q[1];
rz(-0.88417378) q[1];
sx q[1];
rz(-2.2254506) q[1];
x q[2];
rz(2.3362818) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(0.29176816) q[2];
rz(-3.0388888) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(2.5699239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5105195) q[0];
sx q[0];
rz(-1.4820931) q[0];
sx q[0];
rz(-0.53511329) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35927202) q[2];
sx q[2];
rz(-1.1604939) q[2];
sx q[2];
rz(-1.8112195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6211105) q[1];
sx q[1];
rz(-0.74319327) q[1];
sx q[1];
rz(-0.1964257) q[1];
x q[2];
rz(-0.27922697) q[3];
sx q[3];
rz(-0.46776566) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4335094) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(-0.049830534) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(0.18049151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4806992) q[0];
sx q[0];
rz(-1.1850712) q[0];
sx q[0];
rz(-1.4727403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42787243) q[2];
sx q[2];
rz(-2.7395436) q[2];
sx q[2];
rz(-2.7767162) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0879678) q[1];
sx q[1];
rz(-1.8894203) q[1];
sx q[1];
rz(-1.6941403) q[1];
rz(-1.0749712) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(-1.583741) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-0.25156897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0662721) q[0];
sx q[0];
rz(-0.89162725) q[0];
sx q[0];
rz(-0.22596304) q[0];
rz(-pi) q[1];
rz(2.6816363) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(0.83173448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5800036) q[1];
sx q[1];
rz(-1.5434693) q[1];
sx q[1];
rz(2.9461529) q[1];
rz(-2.6336446) q[3];
sx q[3];
rz(-0.93379279) q[3];
sx q[3];
rz(2.7579443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.183737) q[2];
sx q[2];
rz(-2.6001866) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(-1.7306227) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(0.85012287) q[0];
rz(-1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(2.7485671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089038606) q[0];
sx q[0];
rz(-0.35152838) q[0];
sx q[0];
rz(2.0650484) q[0];
x q[1];
rz(1.3156462) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(-1.6659425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0361573) q[1];
sx q[1];
rz(-1.565355) q[1];
sx q[1];
rz(-2.2219707) q[1];
x q[2];
rz(0.29720184) q[3];
sx q[3];
rz(-1.8330049) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3605911) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(-2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-0.0059676776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(-1.1447385) q[0];
rz(-pi) q[1];
rz(-0.75255021) q[2];
sx q[2];
rz(-1.2298905) q[2];
sx q[2];
rz(-0.37170751) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2318864) q[1];
sx q[1];
rz(-1.6344374) q[1];
sx q[1];
rz(2.1257036) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37200125) q[3];
sx q[3];
rz(-1.5270546) q[3];
sx q[3];
rz(-1.7895167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(2.4068508) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(2.1933864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95420102) q[0];
sx q[0];
rz(-1.750995) q[0];
sx q[0];
rz(-0.76822922) q[0];
rz(-pi) q[1];
rz(2.0886705) q[2];
sx q[2];
rz(-1.9621984) q[2];
sx q[2];
rz(-2.1640167) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5612948) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(2.9194174) q[1];
rz(-pi) q[2];
rz(1.4944581) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(-0.59698856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(-0.24027696) q[2];
rz(2.7622973) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(-2.8906524) q[0];
rz(2.8885686) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1044554) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-0.11549581) q[0];
rz(-2.7025931) q[2];
sx q[2];
rz(-0.3294496) q[2];
sx q[2];
rz(-2.0246558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11201227) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(-0.40257247) q[1];
x q[2];
rz(-1.0274067) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(-0.020997626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(2.4582668) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(-1.9911511) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(-2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1922558) q[0];
sx q[0];
rz(-1.0431457) q[0];
sx q[0];
rz(-1.8532182) q[0];
rz(-1.2859225) q[2];
sx q[2];
rz(-1.3956009) q[2];
sx q[2];
rz(-2.3223562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9156956) q[1];
sx q[1];
rz(-0.86663336) q[1];
sx q[1];
rz(-2.8812863) q[1];
rz(0.31302932) q[3];
sx q[3];
rz(-2.1306681) q[3];
sx q[3];
rz(2.5525639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54291723) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1508355) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(-1.5785718) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-2.0832534) q[2];
sx q[2];
rz(-0.28596157) q[2];
sx q[2];
rz(2.7059976) q[2];
rz(1.3068009) q[3];
sx q[3];
rz(-1.055607) q[3];
sx q[3];
rz(-1.9029688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];