OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(5.2360143) q[0];
sx q[0];
rz(9.4935023) q[0];
rz(-1.3955431) q[1];
sx q[1];
rz(-1.5323324) q[1];
sx q[1];
rz(1.2083763) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3291408) q[0];
sx q[0];
rz(-1.4747696) q[0];
sx q[0];
rz(-3.0949233) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5829093) q[2];
sx q[2];
rz(-1.6965869) q[2];
sx q[2];
rz(-2.8319401) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1111869) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(-2.9524132) q[1];
rz(0.73322202) q[3];
sx q[3];
rz(-0.65922046) q[3];
sx q[3];
rz(-2.1472907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.6750083) q[2];
rz(2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(0.74716032) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.117347) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-1.0447964) q[1];
sx q[1];
rz(2.8443964) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575465) q[0];
sx q[0];
rz(-1.4876801) q[0];
sx q[0];
rz(-0.048464171) q[0];
rz(-pi) q[1];
rz(2.6058795) q[2];
sx q[2];
rz(-0.78768724) q[2];
sx q[2];
rz(0.96131575) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.593593) q[1];
sx q[1];
rz(-1.9108692) q[1];
sx q[1];
rz(-1.026457) q[1];
rz(2.7278565) q[3];
sx q[3];
rz(-2.3351151) q[3];
sx q[3];
rz(1.8959351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.979636) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(2.3017853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6162286) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(1.2603659) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8790881) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(1.0002491) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7786583) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(-1.910701) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.967756) q[3];
sx q[3];
rz(-2.1965115) q[3];
sx q[3];
rz(-3.0415149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(-1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.76386219) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(-2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(-0.96558085) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5774028) q[0];
sx q[0];
rz(-1.2947384) q[0];
sx q[0];
rz(2.3645556) q[0];
rz(-pi) q[1];
rz(0.031164073) q[2];
sx q[2];
rz(-2.2989797) q[2];
sx q[2];
rz(2.503501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5406148) q[1];
sx q[1];
rz(-1.317273) q[1];
sx q[1];
rz(-0.032932245) q[1];
x q[2];
rz(0.097316381) q[3];
sx q[3];
rz(-1.684155) q[3];
sx q[3];
rz(1.8196343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(2.4980513) q[3];
sx q[3];
rz(-1.0361592) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(1.2840282) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(-2.0934385) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832344) q[0];
sx q[0];
rz(-2.2336707) q[0];
sx q[0];
rz(-0.62028424) q[0];
rz(-pi) q[1];
rz(2.5571312) q[2];
sx q[2];
rz(-2.2197154) q[2];
sx q[2];
rz(-1.4635758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54742868) q[1];
sx q[1];
rz(-2.8058726) q[1];
sx q[1];
rz(-2.8572542) q[1];
rz(-pi) q[2];
rz(0.89094152) q[3];
sx q[3];
rz(-1.1853301) q[3];
sx q[3];
rz(-2.1894574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2709048) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(0.90083814) q[2];
rz(-2.0488996) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(-1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(-1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(1.2449107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7398864) q[0];
sx q[0];
rz(-2.5283928) q[0];
sx q[0];
rz(0.99341157) q[0];
x q[1];
rz(-2.6558098) q[2];
sx q[2];
rz(-1.8530288) q[2];
sx q[2];
rz(-0.9529875) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60331261) q[1];
sx q[1];
rz(-1.0744175) q[1];
sx q[1];
rz(0.75116317) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8025493) q[3];
sx q[3];
rz(-1.3290805) q[3];
sx q[3];
rz(-1.7237323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9101377) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-0.22496741) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769619) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(1.6784558) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-2.8889012) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1278909) q[0];
sx q[0];
rz(-1.2780006) q[0];
sx q[0];
rz(-3.091759) q[0];
x q[1];
rz(-1.5451317) q[2];
sx q[2];
rz(-0.63458323) q[2];
sx q[2];
rz(-2.6013825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.686324) q[1];
sx q[1];
rz(-0.9372006) q[1];
sx q[1];
rz(-0.92647657) q[1];
x q[2];
rz(0.85315506) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.6097216) q[2];
rz(-1.948471) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-2.0549205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0367592) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(1.6554792) q[0];
rz(-2.8727818) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-2.862646) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19764087) q[0];
sx q[0];
rz(-0.949172) q[0];
sx q[0];
rz(-0.62244121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(2.6627024) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62372359) q[1];
sx q[1];
rz(-1.9451127) q[1];
sx q[1];
rz(-0.75846292) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99854462) q[3];
sx q[3];
rz(-1.6958941) q[3];
sx q[3];
rz(2.3541114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.802861) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(-1.5489244) q[2];
rz(2.0907949) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(-1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.8883702) q[0];
rz(2.5190917) q[1];
sx q[1];
rz(-1.4600735) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86664591) q[0];
sx q[0];
rz(-1.6410315) q[0];
sx q[0];
rz(-2.3932891) q[0];
x q[1];
rz(-2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(2.7149372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5133363) q[1];
sx q[1];
rz(-0.69111052) q[1];
sx q[1];
rz(-1.4555132) q[1];
rz(-pi) q[2];
rz(1.1504437) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(-2.1136485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(2.1257607) q[2];
rz(-2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(-0.36809665) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(0.01424271) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-0.95294398) q[1];
sx q[1];
rz(-1.7600118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4603699) q[0];
sx q[0];
rz(-1.5813706) q[0];
sx q[0];
rz(-1.4282385) q[0];
x q[1];
rz(3.0773452) q[2];
sx q[2];
rz(-2.148743) q[2];
sx q[2];
rz(2.3648175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33243079) q[1];
sx q[1];
rz(-0.5756439) q[1];
sx q[1];
rz(1.3047421) q[1];
rz(-pi) q[2];
rz(0.53346177) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(-1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(2.1515965) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(-0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0902696) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(-1.8016215) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(0.71628911) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(-2.1469231) q[3];
sx q[3];
rz(-2.4352286) q[3];
sx q[3];
rz(2.2027204) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
