OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(2.1265246) q[0];
rz(1.9782344) q[1];
sx q[1];
rz(-0.42127633) q[1];
sx q[1];
rz(1.2383229) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27313156) q[0];
sx q[0];
rz(-0.95306153) q[0];
sx q[0];
rz(2.4988089) q[0];
rz(-1.1467491) q[2];
sx q[2];
rz(-1.2334261) q[2];
sx q[2];
rz(0.81126311) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.048225798) q[1];
sx q[1];
rz(-1.9328572) q[1];
sx q[1];
rz(-2.0680548) q[1];
rz(-pi) q[2];
rz(-1.7389033) q[3];
sx q[3];
rz(-1.013275) q[3];
sx q[3];
rz(1.3489189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6282661) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(-0.70293054) q[2];
rz(-2.2859196) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6913476) q[0];
sx q[0];
rz(-1.0991199) q[0];
sx q[0];
rz(0.74478373) q[0];
rz(-1.5738457) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(2.1994798) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0364931) q[0];
sx q[0];
rz(-0.78723037) q[0];
sx q[0];
rz(-0.28712846) q[0];
x q[1];
rz(-1.4860542) q[2];
sx q[2];
rz(-2.3872445) q[2];
sx q[2];
rz(-0.95065439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84442816) q[1];
sx q[1];
rz(-1.7065863) q[1];
sx q[1];
rz(0.80597767) q[1];
x q[2];
rz(-0.5861939) q[3];
sx q[3];
rz(-1.6209416) q[3];
sx q[3];
rz(1.3231089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-2.0460184) q[2];
rz(1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(1.3862632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0304994) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(1.7362562) q[0];
rz(-1.3966712) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(2.2415846) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25728658) q[0];
sx q[0];
rz(-2.2955756) q[0];
sx q[0];
rz(-1.0777362) q[0];
rz(-pi) q[1];
rz(2.5609803) q[2];
sx q[2];
rz(-1.0913278) q[2];
sx q[2];
rz(-1.978385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72199539) q[1];
sx q[1];
rz(-1.3977786) q[1];
sx q[1];
rz(2.3950775) q[1];
rz(-pi) q[2];
rz(2.4803512) q[3];
sx q[3];
rz(-2.0783278) q[3];
sx q[3];
rz(-2.762261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46093837) q[2];
sx q[2];
rz(-2.9260981) q[2];
sx q[2];
rz(0.94949618) q[2];
rz(-2.5203868) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(-1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42156521) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(-0.048728745) q[0];
rz(0.85982927) q[1];
sx q[1];
rz(-0.33165926) q[1];
sx q[1];
rz(0.31160942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8870114) q[0];
sx q[0];
rz(-1.100228) q[0];
sx q[0];
rz(-2.8356617) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0987058) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(-1.8703415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3032461) q[1];
sx q[1];
rz(-0.7184808) q[1];
sx q[1];
rz(1.9796728) q[1];
x q[2];
rz(0.99589234) q[3];
sx q[3];
rz(-2.6290335) q[3];
sx q[3];
rz(-0.14581524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17778808) q[2];
sx q[2];
rz(-0.21169855) q[2];
sx q[2];
rz(1.4908028) q[2];
rz(-0.35495159) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(1.8420334) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(1.8748913) q[0];
rz(-3.025324) q[1];
sx q[1];
rz(-0.99028844) q[1];
sx q[1];
rz(-1.6273392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86858803) q[0];
sx q[0];
rz(-2.3914861) q[0];
sx q[0];
rz(1.1028034) q[0];
x q[1];
rz(-1.350613) q[2];
sx q[2];
rz(-1.6591958) q[2];
sx q[2];
rz(-2.2184586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4323515) q[1];
sx q[1];
rz(-0.54644692) q[1];
sx q[1];
rz(1.6691471) q[1];
rz(1.5290401) q[3];
sx q[3];
rz(-0.72894086) q[3];
sx q[3];
rz(-1.4545914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9194455) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(-1.9661281) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670052) q[0];
sx q[0];
rz(-0.33127221) q[0];
sx q[0];
rz(-0.40147716) q[0];
rz(1.1270771) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(-2.3289767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82040374) q[0];
sx q[0];
rz(-0.93664737) q[0];
sx q[0];
rz(3.1364596) q[0];
rz(-pi) q[1];
rz(2.4402789) q[2];
sx q[2];
rz(-0.35365401) q[2];
sx q[2];
rz(-0.8071227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8725295) q[1];
sx q[1];
rz(-1.7751964) q[1];
sx q[1];
rz(-2.5356631) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74515588) q[3];
sx q[3];
rz(-2.0260915) q[3];
sx q[3];
rz(-0.32768341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9408985) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(-0.19041348) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.93369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3170526) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(3.0188766) q[0];
rz(-1.1931194) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(2.3741123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240399) q[0];
sx q[0];
rz(-1.5195822) q[0];
sx q[0];
rz(2.2177494) q[0];
rz(-2.9323879) q[2];
sx q[2];
rz(-0.38953094) q[2];
sx q[2];
rz(-0.65702932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40015107) q[1];
sx q[1];
rz(-1.6146683) q[1];
sx q[1];
rz(-0.9167826) q[1];
x q[2];
rz(1.5753059) q[3];
sx q[3];
rz(-1.1645551) q[3];
sx q[3];
rz(2.3227228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7439338) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(1.5896612) q[2];
rz(1.6520366) q[3];
sx q[3];
rz(-1.4068312) q[3];
sx q[3];
rz(1.1154729) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3234696) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-0.37539151) q[0];
rz(-0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323384) q[0];
sx q[0];
rz(-0.94037708) q[0];
sx q[0];
rz(-2.4805043) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.402114) q[2];
sx q[2];
rz(-1.8034435) q[2];
sx q[2];
rz(-2.6017435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.091944253) q[1];
sx q[1];
rz(-2.1990805) q[1];
sx q[1];
rz(-2.6960084) q[1];
x q[2];
rz(2.4861195) q[3];
sx q[3];
rz(-0.87067662) q[3];
sx q[3];
rz(1.9775796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6931307) q[2];
sx q[2];
rz(-1.5631661) q[2];
sx q[2];
rz(0.7473839) q[2];
rz(1.735431) q[3];
sx q[3];
rz(-1.9236671) q[3];
sx q[3];
rz(1.5860484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2324227) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(0.7255834) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.888211) q[1];
sx q[1];
rz(2.5673089) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0030356) q[0];
sx q[0];
rz(-1.7302365) q[0];
sx q[0];
rz(-0.11388679) q[0];
rz(-pi) q[1];
rz(-2.7321283) q[2];
sx q[2];
rz(-2.1019693) q[2];
sx q[2];
rz(0.37497463) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6114823) q[1];
sx q[1];
rz(-2.2951295) q[1];
sx q[1];
rz(-2.9899389) q[1];
rz(-0.58135191) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(2.4331828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1672704) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(-2.4794225) q[2];
rz(-2.3769489) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(1.8173328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8793256) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.7437438) q[0];
rz(-2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(1.3319344) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86517143) q[0];
sx q[0];
rz(-1.4201453) q[0];
sx q[0];
rz(3.1106635) q[0];
rz(-pi) q[1];
rz(0.24253129) q[2];
sx q[2];
rz(-1.9743391) q[2];
sx q[2];
rz(1.5833441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8269013) q[1];
sx q[1];
rz(-1.8002292) q[1];
sx q[1];
rz(1.8073842) q[1];
rz(-pi) q[2];
rz(-2.2374956) q[3];
sx q[3];
rz(-2.3257964) q[3];
sx q[3];
rz(2.286269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8269044) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(0.17987128) q[2];
rz(1.3966857) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(2.9250308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77072813) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(2.6188359) q[1];
sx q[1];
rz(-1.1154543) q[1];
sx q[1];
rz(1.9784068) q[1];
rz(0.63498598) q[2];
sx q[2];
rz(-1.2227092) q[2];
sx q[2];
rz(1.2276445) q[2];
rz(0.64416364) q[3];
sx q[3];
rz(-0.5683036) q[3];
sx q[3];
rz(-0.34556942) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
