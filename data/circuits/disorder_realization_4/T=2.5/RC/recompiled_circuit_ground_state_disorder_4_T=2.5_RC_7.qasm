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
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8684611) q[0];
sx q[0];
rz(-0.95306153) q[0];
sx q[0];
rz(-2.4988089) q[0];
x q[1];
rz(-0.36739393) q[2];
sx q[2];
rz(-1.1720554) q[2];
sx q[2];
rz(-0.90786394) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94436344) q[1];
sx q[1];
rz(-0.60603332) q[1];
sx q[1];
rz(-0.89971772) q[1];
x q[2];
rz(1.7389033) q[3];
sx q[3];
rz(-1.013275) q[3];
sx q[3];
rz(1.7926737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51332659) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(2.4386621) q[2];
rz(-0.85567307) q[3];
sx q[3];
rz(-2.296505) q[3];
sx q[3];
rz(1.8446911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4502451) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(0.74478373) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-0.66190043) q[1];
sx q[1];
rz(0.94211284) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73977913) q[0];
sx q[0];
rz(-1.3688068) q[0];
sx q[0];
rz(-2.3752579) q[0];
rz(-pi) q[1];
rz(1.4860542) q[2];
sx q[2];
rz(-2.3872445) q[2];
sx q[2];
rz(0.95065439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2971645) q[1];
sx q[1];
rz(-1.7065863) q[1];
sx q[1];
rz(-0.80597767) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5553988) q[3];
sx q[3];
rz(-1.5206511) q[3];
sx q[3];
rz(-1.3231089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88910237) q[2];
sx q[2];
rz(-2.1847051) q[2];
sx q[2];
rz(2.0460184) q[2];
rz(1.6628294) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(-1.3862632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-0.11109322) q[0];
sx q[0];
rz(-1.7166623) q[0];
sx q[0];
rz(1.4053364) q[0];
rz(-1.7449215) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(-2.2415846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7177795) q[0];
sx q[0];
rz(-2.2909145) q[0];
sx q[0];
rz(-0.49085842) q[0];
rz(-0.75863691) q[2];
sx q[2];
rz(-2.4066145) q[2];
sx q[2];
rz(-2.1211565) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72199539) q[1];
sx q[1];
rz(-1.3977786) q[1];
sx q[1];
rz(2.3950775) q[1];
x q[2];
rz(-0.7358968) q[3];
sx q[3];
rz(-0.80965878) q[3];
sx q[3];
rz(-2.5084605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6806543) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(0.94949618) q[2];
rz(-0.62120581) q[3];
sx q[3];
rz(-0.92531365) q[3];
sx q[3];
rz(1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200274) q[0];
sx q[0];
rz(-1.255144) q[0];
sx q[0];
rz(0.048728745) q[0];
rz(-2.2817634) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(2.8299832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8870114) q[0];
sx q[0];
rz(-2.0413646) q[0];
sx q[0];
rz(-0.30593095) q[0];
rz(-pi) q[1];
rz(-2.0987058) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(1.2712511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3606122) q[1];
sx q[1];
rz(-0.9223088) q[1];
sx q[1];
rz(-0.33456747) q[1];
rz(2.8446557) q[3];
sx q[3];
rz(-1.9949759) q[3];
sx q[3];
rz(0.78510982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.17778808) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(-1.4908028) q[2];
rz(-2.7866411) q[3];
sx q[3];
rz(-1.7345411) q[3];
sx q[3];
rz(1.8420334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0896924) q[0];
sx q[0];
rz(-0.35744748) q[0];
sx q[0];
rz(-1.8748913) q[0];
rz(0.1162687) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(1.6273392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0564041) q[0];
sx q[0];
rz(-1.2582111) q[0];
sx q[0];
rz(-0.87707918) q[0];
rz(1.1852988) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(2.1182107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4323515) q[1];
sx q[1];
rz(-2.5951457) q[1];
sx q[1];
rz(-1.6691471) q[1];
x q[2];
rz(-0.037260696) q[3];
sx q[3];
rz(-2.2989591) q[3];
sx q[3];
rz(-1.3986349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9194455) q[2];
sx q[2];
rz(-0.47480348) q[2];
sx q[2];
rz(-1.9661281) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(-0.97833943) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745875) q[0];
sx q[0];
rz(-2.8103204) q[0];
sx q[0];
rz(0.40147716) q[0];
rz(-2.0145156) q[1];
sx q[1];
rz(-1.8311484) q[1];
sx q[1];
rz(0.81261596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82040374) q[0];
sx q[0];
rz(-2.2049453) q[0];
sx q[0];
rz(-0.0051330785) q[0];
rz(-pi) q[1];
rz(2.4402789) q[2];
sx q[2];
rz(-0.35365401) q[2];
sx q[2];
rz(-0.8071227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9795831) q[1];
sx q[1];
rz(-2.1623731) q[1];
sx q[1];
rz(1.3237557) q[1];
x q[2];
rz(-2.5162302) q[3];
sx q[3];
rz(-2.2917622) q[3];
sx q[3];
rz(1.4537077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20069417) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(1.5852488) q[2];
rz(-2.9511792) q[3];
sx q[3];
rz(-0.75473458) q[3];
sx q[3];
rz(-1.2079027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82454005) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(1.1931194) q[1];
sx q[1];
rz(-2.2331608) q[1];
sx q[1];
rz(-0.76748031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5106414) q[0];
sx q[0];
rz(-2.2167593) q[0];
sx q[0];
rz(-0.064152282) q[0];
rz(-pi) q[1];
rz(-2.9323879) q[2];
sx q[2];
rz(-2.7520617) q[2];
sx q[2];
rz(0.65702932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1134935) q[1];
sx q[1];
rz(-2.4863247) q[1];
sx q[1];
rz(-1.6428309) q[1];
x q[2];
rz(-0.40624491) q[3];
sx q[3];
rz(-1.5666538) q[3];
sx q[3];
rz(0.7501445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3976589) q[2];
sx q[2];
rz(-3.1199516) q[2];
sx q[2];
rz(-1.5896612) q[2];
rz(-1.4895561) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(2.0261197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3234696) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-2.7662011) q[0];
rz(-2.2581532) q[1];
sx q[1];
rz(-1.5366303) q[1];
sx q[1];
rz(2.4129131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7405072) q[0];
sx q[0];
rz(-2.2621763) q[0];
sx q[0];
rz(-2.2702433) q[0];
rz(-1.402114) q[2];
sx q[2];
rz(-1.3381492) q[2];
sx q[2];
rz(2.6017435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77433005) q[1];
sx q[1];
rz(-0.75241295) q[1];
sx q[1];
rz(-1.03536) q[1];
rz(-pi) q[2];
rz(-2.1971134) q[3];
sx q[3];
rz(-2.2221643) q[3];
sx q[3];
rz(-0.29069549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4484619) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9091699) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(-2.4160093) q[0];
rz(1.8046509) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(0.57428378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1385571) q[0];
sx q[0];
rz(-1.4113562) q[0];
sx q[0];
rz(3.0277059) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7321283) q[2];
sx q[2];
rz(-1.0396233) q[2];
sx q[2];
rz(-0.37497463) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6114823) q[1];
sx q[1];
rz(-0.84646314) q[1];
sx q[1];
rz(-2.9899389) q[1];
x q[2];
rz(-2.1508997) q[3];
sx q[3];
rz(-0.87558666) q[3];
sx q[3];
rz(3.0532072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1672704) q[2];
sx q[2];
rz(-1.7631301) q[2];
sx q[2];
rz(0.66217011) q[2];
rz(-0.76464379) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(1.3242599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.8793256) q[0];
sx q[0];
rz(-2.5536394) q[0];
sx q[0];
rz(-1.3978488) q[0];
rz(2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(1.8096583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66187132) q[0];
sx q[0];
rz(-0.1537696) q[0];
sx q[0];
rz(-1.7717621) q[0];
rz(0.24253129) q[2];
sx q[2];
rz(-1.9743391) q[2];
sx q[2];
rz(-1.5582486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20132682) q[1];
sx q[1];
rz(-1.801071) q[1];
sx q[1];
rz(0.23576945) q[1];
rz(2.5601848) q[3];
sx q[3];
rz(-0.96145844) q[3];
sx q[3];
rz(1.4319624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3146882) q[2];
sx q[2];
rz(-1.1004227) q[2];
sx q[2];
rz(-2.9617214) q[2];
rz(1.3966857) q[3];
sx q[3];
rz(-1.4254009) q[3];
sx q[3];
rz(0.21656187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3708645) q[0];
sx q[0];
rz(-2.6612119) q[0];
sx q[0];
rz(-2.2330855) q[0];
rz(0.52275672) q[1];
sx q[1];
rz(-2.0261384) q[1];
sx q[1];
rz(-1.1631858) q[1];
rz(-2.5066067) q[2];
sx q[2];
rz(-1.2227092) q[2];
sx q[2];
rz(1.2276445) q[2];
rz(0.47209817) q[3];
sx q[3];
rz(-1.2416799) q[3];
sx q[3];
rz(-2.480686) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
