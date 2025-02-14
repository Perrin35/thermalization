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
rz(0.49405721) q[0];
sx q[0];
rz(-0.20209514) q[0];
sx q[0];
rz(-2.8429514) q[0];
rz(2.5477297) q[1];
sx q[1];
rz(-1.9057823) q[1];
sx q[1];
rz(-1.2716582) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32330552) q[0];
sx q[0];
rz(-3.0491017) q[0];
sx q[0];
rz(-2.9748067) q[0];
rz(3.0230672) q[2];
sx q[2];
rz(-1.4473607) q[2];
sx q[2];
rz(0.92230421) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8030111) q[1];
sx q[1];
rz(-2.525248) q[1];
sx q[1];
rz(-2.0508462) q[1];
rz(-pi) q[2];
rz(1.1998953) q[3];
sx q[3];
rz(-1.8270396) q[3];
sx q[3];
rz(-0.56396644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9876081) q[2];
sx q[2];
rz(-2.5488904) q[2];
sx q[2];
rz(0.1565557) q[2];
rz(-0.17078677) q[3];
sx q[3];
rz(-2.4510866) q[3];
sx q[3];
rz(-1.7295711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60468948) q[0];
sx q[0];
rz(-1.5770183) q[0];
sx q[0];
rz(2.1942196) q[0];
rz(2.0715711) q[1];
sx q[1];
rz(-2.1242296) q[1];
sx q[1];
rz(2.1224461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4812093) q[0];
sx q[0];
rz(-2.1041181) q[0];
sx q[0];
rz(-2.9571499) q[0];
x q[1];
rz(1.972946) q[2];
sx q[2];
rz(-2.0216591) q[2];
sx q[2];
rz(0.27489812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85657287) q[1];
sx q[1];
rz(-2.0025221) q[1];
sx q[1];
rz(1.747471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3955973) q[3];
sx q[3];
rz(-2.2713985) q[3];
sx q[3];
rz(1.7563411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77000874) q[2];
sx q[2];
rz(-1.2086478) q[2];
sx q[2];
rz(-2.6303975) q[2];
rz(-0.33358556) q[3];
sx q[3];
rz(-0.89190069) q[3];
sx q[3];
rz(0.79264486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0148575) q[0];
sx q[0];
rz(-0.52913409) q[0];
sx q[0];
rz(1.0134617) q[0];
rz(-0.86499372) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(-1.6046883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3000001) q[0];
sx q[0];
rz(-0.16491297) q[0];
sx q[0];
rz(0.72248904) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75622497) q[2];
sx q[2];
rz(-0.87487223) q[2];
sx q[2];
rz(0.24480259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.930248) q[1];
sx q[1];
rz(-0.67807635) q[1];
sx q[1];
rz(-0.67386595) q[1];
rz(-2.15083) q[3];
sx q[3];
rz(-1.2427075) q[3];
sx q[3];
rz(0.46991959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65261877) q[2];
sx q[2];
rz(-1.3451312) q[2];
sx q[2];
rz(-0.49749231) q[2];
rz(-0.8688212) q[3];
sx q[3];
rz(-2.7610064) q[3];
sx q[3];
rz(3.0697921) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82197613) q[0];
sx q[0];
rz(-3.138534) q[0];
sx q[0];
rz(-2.7724566) q[0];
rz(-0.71040756) q[1];
sx q[1];
rz(-1.0514759) q[1];
sx q[1];
rz(-1.2999473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4669236) q[0];
sx q[0];
rz(-0.82991582) q[0];
sx q[0];
rz(-0.10167565) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7252015) q[2];
sx q[2];
rz(-1.1325628) q[2];
sx q[2];
rz(2.0119502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.877546) q[1];
sx q[1];
rz(-0.49542557) q[1];
sx q[1];
rz(-0.77705125) q[1];
x q[2];
rz(-0.10792144) q[3];
sx q[3];
rz(-2.0156337) q[3];
sx q[3];
rz(0.51714424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15061441) q[2];
sx q[2];
rz(-1.4563478) q[2];
sx q[2];
rz(-0.72188226) q[2];
rz(2.7041096) q[3];
sx q[3];
rz(-2.7719154) q[3];
sx q[3];
rz(-0.30445254) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40780145) q[0];
sx q[0];
rz(-2.3423539) q[0];
sx q[0];
rz(2.9275628) q[0];
rz(2.3389544) q[1];
sx q[1];
rz(-1.9791578) q[1];
sx q[1];
rz(-2.7208557) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20223742) q[0];
sx q[0];
rz(-1.4927683) q[0];
sx q[0];
rz(-1.3952888) q[0];
x q[1];
rz(0.56511648) q[2];
sx q[2];
rz(-1.3656797) q[2];
sx q[2];
rz(2.2895333) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.85462697) q[1];
sx q[1];
rz(-1.4564774) q[1];
sx q[1];
rz(-0.066332093) q[1];
x q[2];
rz(2.2989612) q[3];
sx q[3];
rz(-1.4813322) q[3];
sx q[3];
rz(-2.7842229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3558827) q[2];
sx q[2];
rz(-2.7648401) q[2];
sx q[2];
rz(-0.20998391) q[2];
rz(-2.0436132) q[3];
sx q[3];
rz(-0.72303253) q[3];
sx q[3];
rz(0.22600225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534773) q[0];
sx q[0];
rz(-0.90918875) q[0];
sx q[0];
rz(-0.01532456) q[0];
rz(0.76459908) q[1];
sx q[1];
rz(-1.8761643) q[1];
sx q[1];
rz(-2.9008124) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2368589) q[0];
sx q[0];
rz(-0.52182996) q[0];
sx q[0];
rz(-2.2068992) q[0];
rz(-pi) q[1];
rz(2.639995) q[2];
sx q[2];
rz(-1.7946417) q[2];
sx q[2];
rz(0.9632335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1318553) q[1];
sx q[1];
rz(-2.9564361) q[1];
sx q[1];
rz(-1.2082165) q[1];
rz(-pi) q[2];
rz(0.36846675) q[3];
sx q[3];
rz(-2.0644912) q[3];
sx q[3];
rz(-1.5696978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8814055) q[2];
sx q[2];
rz(-1.0522965) q[2];
sx q[2];
rz(0.50797272) q[2];
rz(2.9925665) q[3];
sx q[3];
rz(-0.38145724) q[3];
sx q[3];
rz(-1.6833359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7070865) q[0];
sx q[0];
rz(-0.91489804) q[0];
sx q[0];
rz(-1.9515422) q[0];
rz(-3.1022364) q[1];
sx q[1];
rz(-2.1085565) q[1];
sx q[1];
rz(-2.8796223) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3888689) q[0];
sx q[0];
rz(-1.5234199) q[0];
sx q[0];
rz(2.8073244) q[0];
rz(-pi) q[1];
rz(-2.8112721) q[2];
sx q[2];
rz(-1.2093423) q[2];
sx q[2];
rz(2.9320331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2319164) q[1];
sx q[1];
rz(-1.2753829) q[1];
sx q[1];
rz(0.40628237) q[1];
rz(-pi) q[2];
rz(2.5983635) q[3];
sx q[3];
rz(-1.3564988) q[3];
sx q[3];
rz(2.8057536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94595853) q[2];
sx q[2];
rz(-2.8527263) q[2];
sx q[2];
rz(-1.757901) q[2];
rz(-1.9687999) q[3];
sx q[3];
rz(-1.5044418) q[3];
sx q[3];
rz(-0.25243944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70343542) q[0];
sx q[0];
rz(-0.96980888) q[0];
sx q[0];
rz(-2.5954212) q[0];
rz(-2.3225885) q[1];
sx q[1];
rz(-1.1823697) q[1];
sx q[1];
rz(-0.46772734) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2692102) q[0];
sx q[0];
rz(-0.2650241) q[0];
sx q[0];
rz(-1.3779968) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1557577) q[2];
sx q[2];
rz(-1.2097683) q[2];
sx q[2];
rz(-2.4280809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7744171) q[1];
sx q[1];
rz(-2.0655334) q[1];
sx q[1];
rz(2.6082188) q[1];
rz(-pi) q[2];
x q[2];
rz(1.813935) q[3];
sx q[3];
rz(-1.7347276) q[3];
sx q[3];
rz(2.1396014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8972682) q[2];
sx q[2];
rz(-1.6017598) q[2];
sx q[2];
rz(-0.77678219) q[2];
rz(1.6952093) q[3];
sx q[3];
rz(-1.8518238) q[3];
sx q[3];
rz(0.18073925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.986213) q[0];
sx q[0];
rz(-1.7387094) q[0];
sx q[0];
rz(-0.078393161) q[0];
rz(-2.8502803) q[1];
sx q[1];
rz(-2.4261609) q[1];
sx q[1];
rz(-1.8142456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9175668) q[0];
sx q[0];
rz(-2.1619658) q[0];
sx q[0];
rz(-2.3252208) q[0];
rz(-0.38703309) q[2];
sx q[2];
rz(-3.0245028) q[2];
sx q[2];
rz(0.69293222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2208682) q[1];
sx q[1];
rz(-2.4358542) q[1];
sx q[1];
rz(-0.42771139) q[1];
x q[2];
rz(-0.45203182) q[3];
sx q[3];
rz(-2.8252183) q[3];
sx q[3];
rz(-1.3133501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3751117) q[2];
sx q[2];
rz(-1.3263308) q[2];
sx q[2];
rz(-0.43409902) q[2];
rz(-2.8050174) q[3];
sx q[3];
rz(-1.429129) q[3];
sx q[3];
rz(-2.7460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340988) q[0];
sx q[0];
rz(-0.33768001) q[0];
sx q[0];
rz(0.39921528) q[0];
rz(-2.4576064) q[1];
sx q[1];
rz(-1.8268879) q[1];
sx q[1];
rz(1.4113873) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1841662) q[0];
sx q[0];
rz(-1.3940576) q[0];
sx q[0];
rz(-2.806163) q[0];
rz(-1.3547276) q[2];
sx q[2];
rz(-2.2377295) q[2];
sx q[2];
rz(1.2264281) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.001907) q[1];
sx q[1];
rz(-1.9314497) q[1];
sx q[1];
rz(2.4649629) q[1];
rz(-2.3606922) q[3];
sx q[3];
rz(-1.077543) q[3];
sx q[3];
rz(-0.79003143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6902249) q[2];
sx q[2];
rz(-1.5294315) q[2];
sx q[2];
rz(-2.8239047) q[2];
rz(2.132527) q[3];
sx q[3];
rz(-1.7694446) q[3];
sx q[3];
rz(2.7528929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379363) q[0];
sx q[0];
rz(-1.6725578) q[0];
sx q[0];
rz(-2.9148711) q[0];
rz(0.33568385) q[1];
sx q[1];
rz(-0.93135584) q[1];
sx q[1];
rz(2.9109536) q[1];
rz(-3.1138218) q[2];
sx q[2];
rz(-0.66412974) q[2];
sx q[2];
rz(0.017505125) q[2];
rz(2.5337467) q[3];
sx q[3];
rz(-1.9030149) q[3];
sx q[3];
rz(-0.72227605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
