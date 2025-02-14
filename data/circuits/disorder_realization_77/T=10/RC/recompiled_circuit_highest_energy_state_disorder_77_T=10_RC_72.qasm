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
rz(-0.51802975) q[0];
sx q[0];
rz(1.1413483) q[0];
sx q[0];
rz(9.4649014) q[0];
rz(-2.8934381) q[1];
sx q[1];
rz(-2.0791972) q[1];
sx q[1];
rz(-3.0575633) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7617159) q[0];
sx q[0];
rz(-1.6296224) q[0];
sx q[0];
rz(0.071278871) q[0];
x q[1];
rz(-2.7438091) q[2];
sx q[2];
rz(-1.4062464) q[2];
sx q[2];
rz(2.0497422) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1511615) q[1];
sx q[1];
rz(-2.6365981) q[1];
sx q[1];
rz(-2.7846365) q[1];
rz(2.5474882) q[3];
sx q[3];
rz(-0.66451529) q[3];
sx q[3];
rz(2.6474109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4485126) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(-2.8903294) q[2];
rz(-1.9692028) q[3];
sx q[3];
rz(-1.1029714) q[3];
sx q[3];
rz(0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6867111) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(-1.7676109) q[0];
rz(-0.33085597) q[1];
sx q[1];
rz(-1.4127981) q[1];
sx q[1];
rz(2.8673598) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5170256) q[0];
sx q[0];
rz(-0.32450482) q[0];
sx q[0];
rz(-0.22402482) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.129446) q[2];
sx q[2];
rz(-0.53330219) q[2];
sx q[2];
rz(-0.27517327) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80503217) q[1];
sx q[1];
rz(-2.7471886) q[1];
sx q[1];
rz(-2.5488362) q[1];
rz(-1.7555439) q[3];
sx q[3];
rz(-0.42506252) q[3];
sx q[3];
rz(-0.63068542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8809044) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(2.0534959) q[2];
rz(-2.0982096) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(2.8823631) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1328717) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(-1.1016499) q[0];
rz(2.3654826) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(-0.64220846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896079) q[0];
sx q[0];
rz(-0.63417681) q[0];
sx q[0];
rz(-0.80192566) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3234336) q[2];
sx q[2];
rz(-2.2752834) q[2];
sx q[2];
rz(0.099964634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1959508) q[1];
sx q[1];
rz(-1.4671051) q[1];
sx q[1];
rz(-1.9983464) q[1];
rz(-pi) q[2];
rz(-2.3989556) q[3];
sx q[3];
rz(-2.3315695) q[3];
sx q[3];
rz(-2.1049316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3099826) q[2];
sx q[2];
rz(-1.1534561) q[2];
sx q[2];
rz(2.2115121) q[2];
rz(0.81405226) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(-1.2306151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26152465) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(-2.6370908) q[0];
rz(0.91744676) q[1];
sx q[1];
rz(-0.73926273) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298131) q[0];
sx q[0];
rz(-1.7195452) q[0];
sx q[0];
rz(2.896033) q[0];
x q[1];
rz(0.31625749) q[2];
sx q[2];
rz(-0.74511792) q[2];
sx q[2];
rz(-0.032501246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.444297) q[1];
sx q[1];
rz(-2.136886) q[1];
sx q[1];
rz(0.56900153) q[1];
rz(-pi) q[2];
rz(-1.3127906) q[3];
sx q[3];
rz(-2.1157221) q[3];
sx q[3];
rz(-1.2598002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.31671277) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(1.971395) q[2];
rz(0.96585387) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-2.476995) q[0];
rz(-0.26847863) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(-2.1427515) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12832175) q[0];
sx q[0];
rz(-1.4185462) q[0];
sx q[0];
rz(0.98441846) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2810318) q[2];
sx q[2];
rz(-2.701607) q[2];
sx q[2];
rz(2.5948465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4376862) q[1];
sx q[1];
rz(-1.4515299) q[1];
sx q[1];
rz(1.4005303) q[1];
rz(-pi) q[2];
rz(-3.1206297) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(1.5251336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1642248) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(-3.0244136) q[2];
rz(0.059727877) q[3];
sx q[3];
rz(-1.2195339) q[3];
sx q[3];
rz(0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(-1.8814948) q[0];
rz(3.0267808) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(3.0462435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95923808) q[0];
sx q[0];
rz(-1.3683912) q[0];
sx q[0];
rz(1.5558045) q[0];
rz(2.9750455) q[2];
sx q[2];
rz(-2.0689365) q[2];
sx q[2];
rz(1.1135676) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9748677) q[1];
sx q[1];
rz(-2.0151867) q[1];
sx q[1];
rz(-0.4687029) q[1];
rz(-pi) q[2];
rz(-0.010706832) q[3];
sx q[3];
rz(-1.043415) q[3];
sx q[3];
rz(-0.37400613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0872588) q[2];
sx q[2];
rz(-0.91780535) q[2];
sx q[2];
rz(-0.66162649) q[2];
rz(2.3719487) q[3];
sx q[3];
rz(-1.6911643) q[3];
sx q[3];
rz(-2.2746287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041644) q[0];
sx q[0];
rz(-2.7935109) q[0];
sx q[0];
rz(-0.75781703) q[0];
rz(3.1163395) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(-1.9662366) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34285082) q[0];
sx q[0];
rz(-1.2161868) q[0];
sx q[0];
rz(0.50574982) q[0];
rz(-pi) q[1];
rz(-0.48904959) q[2];
sx q[2];
rz(-0.95121517) q[2];
sx q[2];
rz(-1.8085768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84056329) q[1];
sx q[1];
rz(-1.1811273) q[1];
sx q[1];
rz(2.3430166) q[1];
rz(-pi) q[2];
rz(-2.712599) q[3];
sx q[3];
rz(-1.5916263) q[3];
sx q[3];
rz(2.5677537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5777099) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(-1.2152524) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(-1.5569713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4828846) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-0.38920745) q[0];
rz(-3.0534741) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22613285) q[0];
sx q[0];
rz(-2.8185039) q[0];
sx q[0];
rz(-2.3996572) q[0];
rz(-pi) q[1];
rz(-1.1168209) q[2];
sx q[2];
rz(-0.96056998) q[2];
sx q[2];
rz(0.35200859) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.033578) q[1];
sx q[1];
rz(-2.696175) q[1];
sx q[1];
rz(3.0207795) q[1];
x q[2];
rz(-2.2810488) q[3];
sx q[3];
rz(-1.8657953) q[3];
sx q[3];
rz(-0.33911946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4072998) q[2];
sx q[2];
rz(-0.6137085) q[2];
sx q[2];
rz(-1.8180004) q[2];
rz(-1.343441) q[3];
sx q[3];
rz(-2.3517793) q[3];
sx q[3];
rz(-2.8048803) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6650498) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(2.6764349) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(2.6506298) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652654) q[0];
sx q[0];
rz(-2.6748228) q[0];
sx q[0];
rz(-0.76961036) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6870306) q[2];
sx q[2];
rz(-2.801932) q[2];
sx q[2];
rz(-0.31955921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1334573) q[1];
sx q[1];
rz(-0.89673364) q[1];
sx q[1];
rz(-2.4070255) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.773953) q[3];
sx q[3];
rz(-0.30045569) q[3];
sx q[3];
rz(-0.57747546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(-2.2779706) q[2];
rz(0.18381271) q[3];
sx q[3];
rz(-0.96537662) q[3];
sx q[3];
rz(2.1553154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154536) q[0];
sx q[0];
rz(-1.621839) q[0];
sx q[0];
rz(-0.62582985) q[0];
rz(2.4848056) q[1];
sx q[1];
rz(-2.1757809) q[1];
sx q[1];
rz(1.9331369) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55126429) q[0];
sx q[0];
rz(-2.9737824) q[0];
sx q[0];
rz(-0.91785927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65932806) q[2];
sx q[2];
rz(-1.8272597) q[2];
sx q[2];
rz(-1.9224482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48599916) q[1];
sx q[1];
rz(-1.5734377) q[1];
sx q[1];
rz(-3.0551617) q[1];
x q[2];
rz(0.27140433) q[3];
sx q[3];
rz(-2.5529666) q[3];
sx q[3];
rz(1.7388104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5586231) q[2];
sx q[2];
rz(-1.8331336) q[2];
sx q[2];
rz(-0.51897955) q[2];
rz(-2.6988622) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(1.7980661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9291572) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(-2.9875372) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(3.0689756) q[2];
sx q[2];
rz(-0.5177064) q[2];
sx q[2];
rz(0.032042423) q[2];
rz(-2.2440425) q[3];
sx q[3];
rz(-1.297045) q[3];
sx q[3];
rz(0.53551077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
