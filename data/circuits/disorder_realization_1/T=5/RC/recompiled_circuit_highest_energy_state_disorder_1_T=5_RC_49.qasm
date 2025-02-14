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
rz(3.0930003) q[0];
sx q[0];
rz(-1.4253923) q[0];
sx q[0];
rz(0.26560321) q[0];
rz(-0.44759294) q[1];
sx q[1];
rz(-1.6375374) q[1];
sx q[1];
rz(-0.12463364) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87089964) q[0];
sx q[0];
rz(-2.2147372) q[0];
sx q[0];
rz(0.55418153) q[0];
x q[1];
rz(1.1963449) q[2];
sx q[2];
rz(-1.6902897) q[2];
sx q[2];
rz(-1.6993285) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9198299) q[1];
sx q[1];
rz(-1.7518834) q[1];
sx q[1];
rz(2.9403375) q[1];
x q[2];
rz(-0.83104317) q[3];
sx q[3];
rz(-0.22749113) q[3];
sx q[3];
rz(0.25615197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0521741) q[2];
sx q[2];
rz(-0.82145059) q[2];
sx q[2];
rz(-2.2649412) q[2];
rz(2.9540201) q[3];
sx q[3];
rz(-2.1552174) q[3];
sx q[3];
rz(-3.0228289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0117884) q[0];
sx q[0];
rz(-0.056948245) q[0];
sx q[0];
rz(2.7563128) q[0];
rz(0.41703364) q[1];
sx q[1];
rz(-2.4200491) q[1];
sx q[1];
rz(-2.1293652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.130439) q[0];
sx q[0];
rz(-2.1277804) q[0];
sx q[0];
rz(0.2574284) q[0];
x q[1];
rz(-0.56372042) q[2];
sx q[2];
rz(-2.4032058) q[2];
sx q[2];
rz(-2.9132879) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31237797) q[1];
sx q[1];
rz(-2.3350403) q[1];
sx q[1];
rz(-0.4017425) q[1];
rz(-2.2506563) q[3];
sx q[3];
rz(-1.9417986) q[3];
sx q[3];
rz(-1.5901417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1483083) q[2];
sx q[2];
rz(-0.40996429) q[2];
sx q[2];
rz(0.5245463) q[2];
rz(-2.2927393) q[3];
sx q[3];
rz(-0.69791228) q[3];
sx q[3];
rz(-0.85936385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9808905) q[0];
sx q[0];
rz(-0.41146678) q[0];
sx q[0];
rz(2.8149862) q[0];
rz(-1.3702565) q[1];
sx q[1];
rz(-1.0561008) q[1];
sx q[1];
rz(-3.1050217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1550138) q[0];
sx q[0];
rz(-0.5875007) q[0];
sx q[0];
rz(-0.67836887) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6447949) q[2];
sx q[2];
rz(-2.3721176) q[2];
sx q[2];
rz(1.066832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0690511) q[1];
sx q[1];
rz(-0.1920192) q[1];
sx q[1];
rz(1.3429696) q[1];
rz(2.1362334) q[3];
sx q[3];
rz(-2.812139) q[3];
sx q[3];
rz(-2.0170596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44564351) q[2];
sx q[2];
rz(-2.8358938) q[2];
sx q[2];
rz(1.9574399) q[2];
rz(-2.6792157) q[3];
sx q[3];
rz(-1.0824883) q[3];
sx q[3];
rz(-2.6760127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9584123) q[0];
sx q[0];
rz(-1.4875702) q[0];
sx q[0];
rz(-2.2271449) q[0];
rz(2.27521) q[1];
sx q[1];
rz(-1.2976846) q[1];
sx q[1];
rz(-2.8241209) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1007967) q[0];
sx q[0];
rz(-1.4616175) q[0];
sx q[0];
rz(-2.7666356) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9725665) q[2];
sx q[2];
rz(-1.9155972) q[2];
sx q[2];
rz(1.930869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9383512) q[1];
sx q[1];
rz(-0.31107357) q[1];
sx q[1];
rz(-2.0955906) q[1];
x q[2];
rz(1.5741979) q[3];
sx q[3];
rz(-0.38013422) q[3];
sx q[3];
rz(1.8148607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.094396599) q[2];
sx q[2];
rz(-2.6698343) q[2];
sx q[2];
rz(-2.6867234) q[2];
rz(-1.1826078) q[3];
sx q[3];
rz(-0.22795658) q[3];
sx q[3];
rz(-0.14127775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-0.32893786) q[0];
sx q[0];
rz(-1.8998572) q[0];
sx q[0];
rz(-0.051359635) q[0];
rz(2.8142169) q[1];
sx q[1];
rz(-2.2739669) q[1];
sx q[1];
rz(0.059965722) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3695196) q[0];
sx q[0];
rz(-0.74893337) q[0];
sx q[0];
rz(-2.0709442) q[0];
rz(-pi) q[1];
rz(-2.4879901) q[2];
sx q[2];
rz(-1.8817668) q[2];
sx q[2];
rz(0.59864908) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9717945) q[1];
sx q[1];
rz(-0.64193875) q[1];
sx q[1];
rz(-1.0331912) q[1];
rz(-0.19503212) q[3];
sx q[3];
rz(-0.86544207) q[3];
sx q[3];
rz(-2.5759199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3943693) q[2];
sx q[2];
rz(-0.9129492) q[2];
sx q[2];
rz(0.24820122) q[2];
rz(-0.40211755) q[3];
sx q[3];
rz(-0.44200236) q[3];
sx q[3];
rz(-0.88002747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1345271) q[0];
sx q[0];
rz(-1.2638673) q[0];
sx q[0];
rz(-1.4136575) q[0];
rz(1.9386579) q[1];
sx q[1];
rz(-0.92034942) q[1];
sx q[1];
rz(-0.10341067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2884379) q[0];
sx q[0];
rz(-1.5010035) q[0];
sx q[0];
rz(1.5091108) q[0];
rz(-pi) q[1];
rz(2.8273583) q[2];
sx q[2];
rz(-1.7467919) q[2];
sx q[2];
rz(2.2361148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0573961) q[1];
sx q[1];
rz(-1.6954719) q[1];
sx q[1];
rz(1.5898819) q[1];
x q[2];
rz(0.2528905) q[3];
sx q[3];
rz(-1.1973945) q[3];
sx q[3];
rz(1.4772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6695413) q[2];
sx q[2];
rz(-2.5504888) q[2];
sx q[2];
rz(-2.9284787) q[2];
rz(-1.537568) q[3];
sx q[3];
rz(-3.0209916) q[3];
sx q[3];
rz(-0.11046256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9977426) q[0];
sx q[0];
rz(-0.85875964) q[0];
sx q[0];
rz(2.6570038) q[0];
rz(-0.46802014) q[1];
sx q[1];
rz(-1.3222398) q[1];
sx q[1];
rz(3.0048634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.384308) q[0];
sx q[0];
rz(-1.1955313) q[0];
sx q[0];
rz(0.25887667) q[0];
rz(-pi) q[1];
rz(-1.9152476) q[2];
sx q[2];
rz(-1.9392804) q[2];
sx q[2];
rz(-2.0257575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.90367442) q[1];
sx q[1];
rz(-0.72924239) q[1];
sx q[1];
rz(-1.7715471) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26094143) q[3];
sx q[3];
rz(-2.709124) q[3];
sx q[3];
rz(-2.3245426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4097269) q[2];
sx q[2];
rz(-3.0146283) q[2];
sx q[2];
rz(-0.48180386) q[2];
rz(-1.9110154) q[3];
sx q[3];
rz(-1.9989719) q[3];
sx q[3];
rz(0.13550152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9442673) q[0];
sx q[0];
rz(-0.33536401) q[0];
sx q[0];
rz(-0.37099567) q[0];
rz(-0.15759298) q[1];
sx q[1];
rz(-1.7540365) q[1];
sx q[1];
rz(-0.88034672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1938095) q[0];
sx q[0];
rz(-2.0447461) q[0];
sx q[0];
rz(3.0633846) q[0];
x q[1];
rz(-1.738638) q[2];
sx q[2];
rz(-0.8693822) q[2];
sx q[2];
rz(2.0244903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.62746) q[1];
sx q[1];
rz(-1.4480906) q[1];
sx q[1];
rz(3.1377327) q[1];
x q[2];
rz(0.28435376) q[3];
sx q[3];
rz(-1.2907249) q[3];
sx q[3];
rz(1.7650616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8627491) q[2];
sx q[2];
rz(-0.66484386) q[2];
sx q[2];
rz(2.1319907) q[2];
rz(-2.3889715) q[3];
sx q[3];
rz(-1.3043343) q[3];
sx q[3];
rz(-2.2032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272334) q[0];
sx q[0];
rz(-0.14150134) q[0];
sx q[0];
rz(-0.046382647) q[0];
rz(1.6498238) q[1];
sx q[1];
rz(-1.3834508) q[1];
sx q[1];
rz(-2.6683064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.804764) q[0];
sx q[0];
rz(-1.7603931) q[0];
sx q[0];
rz(-3.0289344) q[0];
rz(-pi) q[1];
rz(0.32765179) q[2];
sx q[2];
rz(-2.0218414) q[2];
sx q[2];
rz(2.4137677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59236586) q[1];
sx q[1];
rz(-1.8344546) q[1];
sx q[1];
rz(0.63559611) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5978053) q[3];
sx q[3];
rz(-0.031899769) q[3];
sx q[3];
rz(0.75237319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.16372323) q[2];
sx q[2];
rz(-2.551584) q[2];
sx q[2];
rz(-3.1150277) q[2];
rz(2.6350392) q[3];
sx q[3];
rz(-2.3285464) q[3];
sx q[3];
rz(2.626239) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16515054) q[0];
sx q[0];
rz(-0.9557752) q[0];
sx q[0];
rz(0.15560786) q[0];
rz(0.43442976) q[1];
sx q[1];
rz(-0.54672086) q[1];
sx q[1];
rz(2.853552) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1032938) q[0];
sx q[0];
rz(-1.0464051) q[0];
sx q[0];
rz(0.68322261) q[0];
rz(-1.4819809) q[2];
sx q[2];
rz(-0.15379158) q[2];
sx q[2];
rz(1.3336934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.23032863) q[1];
sx q[1];
rz(-0.9254749) q[1];
sx q[1];
rz(-0.55379587) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41573306) q[3];
sx q[3];
rz(-0.98728335) q[3];
sx q[3];
rz(0.34900591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7274999) q[2];
sx q[2];
rz(-0.45656559) q[2];
sx q[2];
rz(-0.045819316) q[2];
rz(3.0231061) q[3];
sx q[3];
rz(-0.89209569) q[3];
sx q[3];
rz(2.4011325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.4988149) q[0];
sx q[0];
rz(-1.5121664) q[0];
sx q[0];
rz(1.2335516) q[0];
rz(1.9998101) q[1];
sx q[1];
rz(-2.4363166) q[1];
sx q[1];
rz(1.8038149) q[1];
rz(2.4579688) q[2];
sx q[2];
rz(-1.0928921) q[2];
sx q[2];
rz(3.0047807) q[2];
rz(1.6330887) q[3];
sx q[3];
rz(-1.6332165) q[3];
sx q[3];
rz(-1.6467057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
