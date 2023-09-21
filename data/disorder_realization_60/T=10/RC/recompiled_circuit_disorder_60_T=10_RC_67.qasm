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
rz(-1.047171) q[0];
sx q[0];
rz(0.068724364) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8954593) q[0];
sx q[0];
rz(-1.6172505) q[0];
sx q[0];
rz(-1.6669271) q[0];
rz(-pi) q[1];
rz(-3/(10*pi)) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(-0.40590826) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1111869) q[1];
sx q[1];
rz(-2.5136247) q[1];
sx q[1];
rz(-2.9524132) q[1];
rz(-2.049202) q[3];
sx q[3];
rz(-2.0432825) q[3];
sx q[3];
rz(2.9977968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(-2.4438434) q[3];
sx q[3];
rz(-1.1013228) q[3];
sx q[3];
rz(2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.117347) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(-1.1741937) q[0];
rz(-2.970447) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78272351) q[0];
sx q[0];
rz(-1.5224996) q[0];
sx q[0];
rz(1.6540098) q[0];
rz(-2.0446288) q[2];
sx q[2];
rz(-2.2261438) q[2];
sx q[2];
rz(-1.4807793) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.593593) q[1];
sx q[1];
rz(-1.2307234) q[1];
sx q[1];
rz(-1.026457) q[1];
x q[2];
rz(-1.1737212) q[3];
sx q[3];
rz(-0.84871549) q[3];
sx q[3];
rz(1.811036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.979636) q[2];
sx q[2];
rz(-1.9571783) q[2];
sx q[2];
rz(2.5276108) q[2];
rz(-0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-0.30763787) q[0];
rz(2.3930507) q[1];
sx q[1];
rz(-0.33154878) q[1];
sx q[1];
rz(0.83980733) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15159431) q[0];
sx q[0];
rz(-2.6322106) q[0];
sx q[0];
rz(2.529782) q[0];
rz(1.4762127) q[2];
sx q[2];
rz(-1.2328096) q[2];
sx q[2];
rz(-2.4199977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3629344) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(1.2308916) q[1];
x q[2];
rz(1.8057459) q[3];
sx q[3];
rz(-2.4953105) q[3];
sx q[3];
rz(2.7502053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86747375) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(1.8236558) q[2];
rz(1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(1.6320451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777305) q[0];
sx q[0];
rz(-1.7611935) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(-0.96558085) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8646116) q[0];
sx q[0];
rz(-0.81482139) q[0];
sx q[0];
rz(2.7576202) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1104286) q[2];
sx q[2];
rz(-2.2989797) q[2];
sx q[2];
rz(2.503501) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6009778) q[1];
sx q[1];
rz(-1.8243196) q[1];
sx q[1];
rz(-0.032932245) q[1];
rz(-pi) q[2];
rz(2.2772917) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(-2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(1.003456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716361) q[0];
sx q[0];
rz(-1.5820553) q[0];
sx q[0];
rz(-1.2840282) q[0];
rz(-2.8517826) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5149887) q[0];
sx q[0];
rz(-2.0467313) q[0];
sx q[0];
rz(-0.80608741) q[0];
rz(-pi) q[1];
rz(2.3088147) q[2];
sx q[2];
rz(-2.0259498) q[2];
sx q[2];
rz(2.8684794) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84761274) q[1];
sx q[1];
rz(-1.8925397) q[1];
sx q[1];
rz(-1.4732248) q[1];
x q[2];
rz(-2.2506511) q[3];
sx q[3];
rz(-1.9562625) q[3];
sx q[3];
rz(-0.95213529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(0.90083814) q[2];
rz(-1.0926931) q[3];
sx q[3];
rz(-1.0034424) q[3];
sx q[3];
rz(1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95935217) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-2.545488) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(-1.2449107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67958528) q[0];
sx q[0];
rz(-1.2512659) q[0];
sx q[0];
rz(-2.1035478) q[0];
rz(-pi) q[1];
rz(-1.8876569) q[2];
sx q[2];
rz(-2.0358026) q[2];
sx q[2];
rz(-2.3777865) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7024755) q[1];
sx q[1];
rz(-2.2687952) q[1];
sx q[1];
rz(2.4707787) q[1];
rz(-pi) q[2];
rz(2.8935029) q[3];
sx q[3];
rz(-1.7956942) q[3];
sx q[3];
rz(-0.096506462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.231455) q[2];
sx q[2];
rz(-2.4427876) q[2];
sx q[2];
rz(2.9166252) q[2];
rz(0.088430017) q[3];
sx q[3];
rz(-1.6902573) q[3];
sx q[3];
rz(0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46463075) q[0];
sx q[0];
rz(-1.9043652) q[0];
sx q[0];
rz(-0.27994573) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(-0.25269145) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7130816) q[0];
sx q[0];
rz(-1.6185074) q[0];
sx q[0];
rz(-1.8639355) q[0];
x q[1];
rz(-1.5451317) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(2.6013825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6934769) q[1];
sx q[1];
rz(-2.2709393) q[1];
sx q[1];
rz(-2.4561988) q[1];
rz(-pi) q[2];
rz(2.3971862) q[3];
sx q[3];
rz(-0.91074569) q[3];
sx q[3];
rz(-0.50124121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10483345) q[0];
sx q[0];
rz(-1.3830673) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(0.2688109) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-2.862646) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689965) q[0];
sx q[0];
rz(-2.0645752) q[0];
sx q[0];
rz(-0.84817024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8770304) q[2];
sx q[2];
rz(-2.6191204) q[2];
sx q[2];
rz(0.82447169) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62372359) q[1];
sx q[1];
rz(-1.9451127) q[1];
sx q[1];
rz(-0.75846292) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9931077) q[3];
sx q[3];
rz(-1.0035702) q[3];
sx q[3];
rz(2.4384769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.4638476) q[2];
sx q[2];
rz(1.5926682) q[2];
rz(-2.0907949) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.9416434) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46362296) q[0];
sx q[0];
rz(-2.2080053) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3619986) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(3.0385618) q[0];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(-2.7149372) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4790926) q[1];
sx q[1];
rz(-2.2564285) q[1];
sx q[1];
rz(-0.094866026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9911489) q[3];
sx q[3];
rz(-0.97878362) q[3];
sx q[3];
rz(-1.0279442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(-3.1273499) q[0];
rz(0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.3815809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89109126) q[0];
sx q[0];
rz(-1.4282465) q[0];
sx q[0];
rz(-0.010682627) q[0];
rz(3.0773452) q[2];
sx q[2];
rz(-2.148743) q[2];
sx q[2];
rz(2.3648175) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.018316293) q[1];
sx q[1];
rz(-2.1237719) q[1];
sx q[1];
rz(-0.1690013) q[1];
rz(0.53346177) q[3];
sx q[3];
rz(-1.3812314) q[3];
sx q[3];
rz(1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(2.1515965) q[2];
rz(2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.0902696) q[0];
sx q[0];
rz(-0.6482424) q[0];
sx q[0];
rz(2.0363664) q[0];
rz(-1.3399711) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(2.6132265) q[2];
sx q[2];
rz(-2.3522204) q[2];
sx q[2];
rz(1.5632202) q[2];
rz(0.99466956) q[3];
sx q[3];
rz(-2.4352286) q[3];
sx q[3];
rz(2.2027204) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
