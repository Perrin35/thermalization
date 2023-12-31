OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(-1.6337448) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36131418) q[0];
sx q[0];
rz(-1.2957934) q[0];
sx q[0];
rz(1.6367903) q[0];
x q[1];
rz(1.3747896) q[2];
sx q[2];
rz(-2.4110846) q[2];
sx q[2];
rz(2.1240049) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7826739) q[1];
sx q[1];
rz(-2.3896857) q[1];
sx q[1];
rz(2.9208675) q[1];
rz(1.4673759) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-2.1372674) q[0];
rz(1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(2.1388334) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28264499) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(2.100201) q[0];
rz(-2.9687256) q[2];
sx q[2];
rz(-1.3860518) q[2];
sx q[2];
rz(2.6181521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7748389) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(-1.8444091) q[1];
rz(-pi) q[2];
rz(-0.49591222) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(1.0752614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-0.77899581) q[0];
rz(-0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(0.88358203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040366216) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(-2.7357487) q[0];
x q[1];
rz(-0.94647879) q[2];
sx q[2];
rz(-0.4779856) q[2];
sx q[2];
rz(-1.1849272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(1.6509389) q[1];
rz(0.15150841) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(-0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7097968) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(-1.1071831) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54116762) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(3.0279972) q[0];
rz(-pi) q[1];
rz(-2.8407211) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(3.0645264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(-0.99460852) q[1];
rz(-pi) q[2];
rz(0.21869603) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(-2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52175534) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0147858) q[0];
sx q[0];
rz(-0.94313398) q[0];
sx q[0];
rz(3.0833551) q[0];
x q[1];
rz(-3.0643164) q[2];
sx q[2];
rz(-0.52646598) q[2];
sx q[2];
rz(-2.4174945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.535714) q[1];
sx q[1];
rz(-1.5389581) q[1];
sx q[1];
rz(-1.72623) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3239273) q[3];
sx q[3];
rz(-0.92901232) q[3];
sx q[3];
rz(-1.6587917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0267462) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-0.39658305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30124861) q[0];
sx q[0];
rz(-2.6123971) q[0];
sx q[0];
rz(-2.089414) q[0];
rz(-pi) q[1];
rz(2.5846892) q[2];
sx q[2];
rz(-0.54585005) q[2];
sx q[2];
rz(2.5715373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0551758) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(-0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(0.82715183) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267245) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(0.29863775) q[0];
rz(-pi) q[1];
rz(2.9992505) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(3.1359283) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78325242) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.7794442) q[1];
rz(-1.9184291) q[3];
sx q[3];
rz(-1.6602483) q[3];
sx q[3];
rz(-1.7355433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(0.032756068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2643124) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(2.2420922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1342437) q[2];
sx q[2];
rz(-1.9512366) q[2];
sx q[2];
rz(-3.0055339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0667463) q[1];
sx q[1];
rz(-0.42145887) q[1];
sx q[1];
rz(-1.4261817) q[1];
rz(-pi) q[2];
rz(1.1912187) q[3];
sx q[3];
rz(-1.6536342) q[3];
sx q[3];
rz(-2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63472414) q[0];
sx q[0];
rz(-1.5171432) q[0];
sx q[0];
rz(-0.13454484) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59818563) q[2];
sx q[2];
rz(-1.7028114) q[2];
sx q[2];
rz(2.8796632) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-2.328518) q[1];
rz(-2.5599307) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(1.7285085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0316524) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1372676) q[0];
sx q[0];
rz(-0.30170545) q[0];
sx q[0];
rz(0.76122491) q[0];
rz(-pi) q[1];
rz(-1.7190476) q[2];
sx q[2];
rz(-1.3607927) q[2];
sx q[2];
rz(1.6809747) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88565247) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(-3.0401858) q[1];
x q[2];
rz(-2.3365799) q[3];
sx q[3];
rz(-0.22634889) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(2.9700206) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(-2.0486352) q[3];
sx q[3];
rz(-2.408705) q[3];
sx q[3];
rz(-0.68791289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
