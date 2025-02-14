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
rz(-3.0012335) q[0];
sx q[0];
rz(-1.4599414) q[0];
sx q[0];
rz(2.2886544) q[0];
rz(-1.1879022) q[1];
sx q[1];
rz(-0.24628425) q[1];
sx q[1];
rz(-2.8077717) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4720195) q[0];
sx q[0];
rz(-1.8209576) q[0];
sx q[0];
rz(-1.0521558) q[0];
rz(-1.0707955) q[2];
sx q[2];
rz(-0.99783932) q[2];
sx q[2];
rz(-0.9129325) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94813628) q[1];
sx q[1];
rz(-1.4206593) q[1];
sx q[1];
rz(2.9562922) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4859441) q[3];
sx q[3];
rz(-1.4699664) q[3];
sx q[3];
rz(-2.3840897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9217801) q[2];
sx q[2];
rz(-0.87679902) q[2];
sx q[2];
rz(-2.5019257) q[2];
rz(-2.2031247) q[3];
sx q[3];
rz(-0.42249051) q[3];
sx q[3];
rz(-1.2708906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6476145) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(-1.6139503) q[0];
rz(-0.24213067) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(0.72227824) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764389) q[0];
sx q[0];
rz(-1.1761888) q[0];
sx q[0];
rz(-2.7662686) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8794888) q[2];
sx q[2];
rz(-2.3293709) q[2];
sx q[2];
rz(1.1103528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3434365) q[1];
sx q[1];
rz(-1.0381471) q[1];
sx q[1];
rz(0.58782593) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8037968) q[3];
sx q[3];
rz(-2.5058441) q[3];
sx q[3];
rz(-1.1274757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.051141288) q[2];
sx q[2];
rz(-2.800056) q[2];
sx q[2];
rz(-0.086645834) q[2];
rz(0.58049774) q[3];
sx q[3];
rz(-1.9167506) q[3];
sx q[3];
rz(2.1897924) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314303) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(-1.6828368) q[0];
rz(0.1156062) q[1];
sx q[1];
rz(-1.1980201) q[1];
sx q[1];
rz(-2.9748532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939759) q[0];
sx q[0];
rz(-2.1361094) q[0];
sx q[0];
rz(0.24576776) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9405792) q[2];
sx q[2];
rz(-1.175011) q[2];
sx q[2];
rz(-1.41768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84126842) q[1];
sx q[1];
rz(-2.1144697) q[1];
sx q[1];
rz(0.44692301) q[1];
rz(-pi) q[2];
rz(2.3963967) q[3];
sx q[3];
rz(-1.6752073) q[3];
sx q[3];
rz(-0.94693434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.705767) q[2];
sx q[2];
rz(-2.9220118) q[2];
sx q[2];
rz(-2.0078008) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.2727126) q[3];
sx q[3];
rz(0.72499609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0024006) q[0];
sx q[0];
rz(-0.63183689) q[0];
sx q[0];
rz(-2.8644526) q[0];
rz(-0.78284043) q[1];
sx q[1];
rz(-1.6354086) q[1];
sx q[1];
rz(-0.35071075) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4970873) q[0];
sx q[0];
rz(-0.86928029) q[0];
sx q[0];
rz(0.10414609) q[0];
rz(0.088557505) q[2];
sx q[2];
rz(-2.8148201) q[2];
sx q[2];
rz(1.8984924) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33138946) q[1];
sx q[1];
rz(-2.1955829) q[1];
sx q[1];
rz(-0.58831711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043204149) q[3];
sx q[3];
rz(-1.6790048) q[3];
sx q[3];
rz(-1.6668591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0413401) q[2];
sx q[2];
rz(-1.2790054) q[2];
sx q[2];
rz(-1.1502385) q[2];
rz(2.9669115) q[3];
sx q[3];
rz(-2.8476069) q[3];
sx q[3];
rz(0.090350769) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028285099) q[0];
sx q[0];
rz(-0.57357967) q[0];
sx q[0];
rz(1.9453402) q[0];
rz(2.664227) q[1];
sx q[1];
rz(-0.82672516) q[1];
sx q[1];
rz(0.54944077) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73348896) q[0];
sx q[0];
rz(-0.86015742) q[0];
sx q[0];
rz(-2.3657794) q[0];
x q[1];
rz(-0.35713335) q[2];
sx q[2];
rz(-2.5741842) q[2];
sx q[2];
rz(2.0878938) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3122299) q[1];
sx q[1];
rz(-2.7155502) q[1];
sx q[1];
rz(-1.2755574) q[1];
rz(-2.7189452) q[3];
sx q[3];
rz(-1.1492582) q[3];
sx q[3];
rz(-2.1545422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4598733) q[2];
sx q[2];
rz(-0.37178603) q[2];
sx q[2];
rz(-2.8968774) q[2];
rz(1.5185482) q[3];
sx q[3];
rz(-1.1145498) q[3];
sx q[3];
rz(-3.0486619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.86843425) q[0];
sx q[0];
rz(-0.99240977) q[0];
sx q[0];
rz(0.42700818) q[0];
rz(1.931841) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(0.0078113656) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46011074) q[0];
sx q[0];
rz(-0.78435271) q[0];
sx q[0];
rz(-1.2624692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1179256) q[2];
sx q[2];
rz(-1.5627699) q[2];
sx q[2];
rz(-3.0723177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4512109) q[1];
sx q[1];
rz(-0.65245095) q[1];
sx q[1];
rz(3.058372) q[1];
rz(0.2711556) q[3];
sx q[3];
rz(-0.89224766) q[3];
sx q[3];
rz(-3.1377047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4012287) q[2];
sx q[2];
rz(-0.87330356) q[2];
sx q[2];
rz(-2.002423) q[2];
rz(-0.27979699) q[3];
sx q[3];
rz(-1.8179025) q[3];
sx q[3];
rz(1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4918268) q[0];
sx q[0];
rz(-0.38182807) q[0];
sx q[0];
rz(2.5873798) q[0];
rz(-0.062049374) q[1];
sx q[1];
rz(-2.4833312) q[1];
sx q[1];
rz(2.2668692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85200441) q[0];
sx q[0];
rz(-1.6920493) q[0];
sx q[0];
rz(-2.3512588) q[0];
rz(-pi) q[1];
rz(-0.0055829835) q[2];
sx q[2];
rz(-1.5632331) q[2];
sx q[2];
rz(0.39952229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8043912) q[1];
sx q[1];
rz(-2.534165) q[1];
sx q[1];
rz(1.1056221) q[1];
x q[2];
rz(1.0933541) q[3];
sx q[3];
rz(-2.4191609) q[3];
sx q[3];
rz(0.6078161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7265861) q[2];
sx q[2];
rz(-1.0174624) q[2];
sx q[2];
rz(0.93878186) q[2];
rz(-2.4472661) q[3];
sx q[3];
rz(-2.4735579) q[3];
sx q[3];
rz(-0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.071455) q[0];
sx q[0];
rz(-2.5818765) q[0];
sx q[0];
rz(0.30858421) q[0];
rz(1.6584819) q[1];
sx q[1];
rz(-1.3456656) q[1];
sx q[1];
rz(2.9187091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6312941) q[0];
sx q[0];
rz(-1.1985072) q[0];
sx q[0];
rz(0.75081236) q[0];
x q[1];
rz(2.9176788) q[2];
sx q[2];
rz(-1.4707047) q[2];
sx q[2];
rz(1.3097641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1403583) q[1];
sx q[1];
rz(-1.3257294) q[1];
sx q[1];
rz(0.040018572) q[1];
rz(-pi) q[2];
rz(-1.0142266) q[3];
sx q[3];
rz(-0.96650079) q[3];
sx q[3];
rz(-0.44318553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8181204) q[2];
sx q[2];
rz(-1.0264531) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.1431665) q[3];
sx q[3];
rz(-1.4799456) q[3];
sx q[3];
rz(-2.0460879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185862) q[0];
sx q[0];
rz(-1.9976595) q[0];
sx q[0];
rz(-3.123172) q[0];
rz(0.1704692) q[1];
sx q[1];
rz(-1.2455995) q[1];
sx q[1];
rz(2.9439994) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9182943) q[0];
sx q[0];
rz(-0.16897783) q[0];
sx q[0];
rz(-1.6264982) q[0];
x q[1];
rz(1.9371447) q[2];
sx q[2];
rz(-1.6197259) q[2];
sx q[2];
rz(1.6530703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8325628) q[1];
sx q[1];
rz(-1.4352928) q[1];
sx q[1];
rz(-1.871528) q[1];
rz(-pi) q[2];
rz(-0.39579795) q[3];
sx q[3];
rz(-1.6230109) q[3];
sx q[3];
rz(2.8767916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0240872) q[2];
sx q[2];
rz(-2.0427637) q[2];
sx q[2];
rz(2.7433024) q[2];
rz(2.6309218) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(-0.85810703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9061822) q[0];
sx q[0];
rz(-0.76950961) q[0];
sx q[0];
rz(2.9421575) q[0];
rz(2.8681352) q[1];
sx q[1];
rz(-0.43379915) q[1];
sx q[1];
rz(1.010703) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.511925) q[0];
sx q[0];
rz(-0.92610794) q[0];
sx q[0];
rz(-1.900338) q[0];
x q[1];
rz(-1.5908786) q[2];
sx q[2];
rz(-1.3354567) q[2];
sx q[2];
rz(-2.2724336) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5291374) q[1];
sx q[1];
rz(-0.83548629) q[1];
sx q[1];
rz(1.6409954) q[1];
rz(-pi) q[2];
rz(-2.3673986) q[3];
sx q[3];
rz(-2.7280118) q[3];
sx q[3];
rz(-0.78998297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79591862) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(0.79552135) q[2];
rz(-0.37122053) q[3];
sx q[3];
rz(-1.3859387) q[3];
sx q[3];
rz(2.8871239) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026074792) q[0];
sx q[0];
rz(-1.4160897) q[0];
sx q[0];
rz(-1.8427451) q[0];
rz(-0.55116354) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(-1.3462832) q[2];
sx q[2];
rz(-2.2835352) q[2];
sx q[2];
rz(2.78751) q[2];
rz(1.9528648) q[3];
sx q[3];
rz(-1.794906) q[3];
sx q[3];
rz(1.256663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
