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
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9809496) q[0];
sx q[0];
rz(-2.0619218) q[0];
sx q[0];
rz(-2.8777697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.067692368) q[2];
sx q[2];
rz(-0.86172416) q[2];
sx q[2];
rz(-0.46082218) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(-2.5488528) q[1];
rz(-pi) q[2];
rz(-2.4077971) q[3];
sx q[3];
rz(-1.9473837) q[3];
sx q[3];
rz(0.22613444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4094231) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(-2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10193292) q[0];
sx q[0];
rz(-1.9358557) q[0];
sx q[0];
rz(1.2714766) q[0];
x q[1];
rz(-0.14658908) q[2];
sx q[2];
rz(-0.97448889) q[2];
sx q[2];
rz(1.5335611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12304141) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(0.80177387) q[1];
x q[2];
rz(2.3362818) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-2.2413072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.7385959) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(2.8318751) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6310731) q[0];
sx q[0];
rz(-1.6594995) q[0];
sx q[0];
rz(-0.53511329) q[0];
rz(-pi) q[1];
rz(-1.13582) q[2];
sx q[2];
rz(-1.8990714) q[2];
sx q[2];
rz(-3.0498743) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6211105) q[1];
sx q[1];
rz(-2.3983994) q[1];
sx q[1];
rz(-0.1964257) q[1];
rz(0.45205558) q[3];
sx q[3];
rz(-1.4462024) q[3];
sx q[3];
rz(0.053754036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9812575) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4335094) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(3.0917621) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-0.18049151) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1946963) q[0];
sx q[0];
rz(-1.4799656) q[0];
sx q[0];
rz(-2.7541861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7137202) q[2];
sx q[2];
rz(-2.7395436) q[2];
sx q[2];
rz(-0.36487647) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6199477) q[1];
sx q[1];
rz(-1.6879028) q[1];
sx q[1];
rz(0.32089969) q[1];
x q[2];
rz(-2.0666215) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(-1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(-1.583741) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(0.98651648) q[0];
rz(-1.9793234) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7154327) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(1.3000814) q[0];
rz(-1.4839843) q[2];
sx q[2];
rz(-1.112339) q[2];
sx q[2];
rz(-2.364033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1269762) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(-1.5986534) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1523347) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(1.1359515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(-0.63123909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(-1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(2.7485671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0133936) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(1.8830995) q[0];
rz(1.3156462) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(-1.6659425) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45822083) q[1];
sx q[1];
rz(-0.65119377) q[1];
sx q[1];
rz(1.5797735) q[1];
rz(-pi) q[2];
rz(1.2971446) q[3];
sx q[3];
rz(-1.8575462) q[3];
sx q[3];
rz(-1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3605911) q[2];
sx q[2];
rz(-1.1258619) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(2.1940103) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.9472286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(-0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3403444) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(1.1447385) q[0];
x q[1];
rz(1.1184095) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(-1.6391022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2318864) q[1];
sx q[1];
rz(-1.6344374) q[1];
sx q[1];
rz(2.1257036) q[1];
x q[2];
rz(-1.5238477) q[3];
sx q[3];
rz(-1.9424244) q[3];
sx q[3];
rz(-0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(2.9515284) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-0.53623143) q[0];
rz(-2.2082632) q[1];
sx q[1];
rz(-1.5737165) q[1];
sx q[1];
rz(0.94820625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78806879) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(1.3226932) q[0];
rz(-2.6981632) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(-2.7623917) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5612948) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(2.9194174) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7847399) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(0.94716351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.5189974) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(0.25094029) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-2.7889263) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0371373) q[0];
sx q[0];
rz(-1.6933352) q[0];
sx q[0];
rz(-0.11549581) q[0];
x q[1];
rz(0.30013957) q[2];
sx q[2];
rz(-1.7087414) q[2];
sx q[2];
rz(-2.2696242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36665146) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(1.304438) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6412233) q[3];
sx q[3];
rz(-0.89958588) q[3];
sx q[3];
rz(-0.74218111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(2.4582668) q[2];
rz(-1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7651354) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(-2.4328649) q[1];
sx q[1];
rz(-0.47660247) q[1];
sx q[1];
rz(-0.93651071) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6652128) q[0];
sx q[0];
rz(-1.8139651) q[0];
sx q[0];
rz(0.54540821) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1328743) q[2];
sx q[2];
rz(-2.808411) q[2];
sx q[2];
rz(2.9269232) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5156538) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(-2.291912) q[1];
rz(-pi) q[2];
rz(-2.0274721) q[3];
sx q[3];
rz(-2.5084087) q[3];
sx q[3];
rz(3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(-2.9144918) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(0.14317748) q[2];
sx q[2];
rz(-1.3224052) q[2];
sx q[2];
rz(0.094766141) q[2];
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
