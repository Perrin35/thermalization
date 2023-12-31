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
rz(-1.9332164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24613334) q[0];
sx q[0];
rz(-1.5243422) q[0];
sx q[0];
rz(-1.6669271) q[0];
rz(1.5586833) q[2];
sx q[2];
rz(-1.4450057) q[2];
sx q[2];
rz(-2.8319401) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0304058) q[1];
sx q[1];
rz(-0.62796794) q[1];
sx q[1];
rz(-0.18917947) q[1];
rz(-pi) q[2];
rz(2.4083706) q[3];
sx q[3];
rz(-2.4823722) q[3];
sx q[3];
rz(0.99430195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3258813) q[2];
sx q[2];
rz(-1.4006961) q[2];
sx q[2];
rz(-1.6750083) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(-0.74716032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0242457) q[0];
sx q[0];
rz(-1.1276561) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(-0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(2.8443964) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78404616) q[0];
sx q[0];
rz(-1.6539126) q[0];
sx q[0];
rz(-3.0931285) q[0];
rz(-2.4291123) q[2];
sx q[2];
rz(-1.2006294) q[2];
sx q[2];
rz(0.21288255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(0.39217197) q[1];
rz(-2.7278565) q[3];
sx q[3];
rz(-2.3351151) q[3];
sx q[3];
rz(1.2456576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.979636) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(0.61398181) q[2];
rz(2.2654514) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(2.5045625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0456332) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(0.30763787) q[0];
rz(0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-2.3017853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1727985) q[0];
sx q[0];
rz(-1.286924) q[0];
sx q[0];
rz(2.7127405) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26250458) q[2];
sx q[2];
rz(-2.7911107) q[2];
sx q[2];
rz(-2.1413435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0963904) q[1];
sx q[1];
rz(-0.55711105) q[1];
sx q[1];
rz(0.60369173) q[1];
rz(-pi) q[2];
rz(2.2037376) q[3];
sx q[3];
rz(-1.4301392) q[3];
sx q[3];
rz(1.7733639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2741189) q[2];
sx q[2];
rz(-1.7904736) q[2];
sx q[2];
rz(-1.3179368) q[2];
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
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777305) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(2.6960301) q[0];
rz(2.5207649) q[1];
sx q[1];
rz(-1.2460243) q[1];
sx q[1];
rz(-2.1760118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25528204) q[0];
sx q[0];
rz(-2.3112486) q[0];
sx q[0];
rz(-1.9489954) q[0];
rz(-pi) q[1];
rz(0.84237174) q[2];
sx q[2];
rz(-1.594055) q[2];
sx q[2];
rz(-2.1881441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4703776) q[1];
sx q[1];
rz(-0.2556076) q[1];
sx q[1];
rz(-1.6971991) q[1];
x q[2];
rz(-0.86430092) q[3];
sx q[3];
rz(-2.9923277) q[3];
sx q[3];
rz(-2.5316558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3199557) q[2];
sx q[2];
rz(-2.379202) q[2];
sx q[2];
rz(-2.2122673) q[2];
rz(0.6435414) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(1.003456) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699566) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(-1.8575645) q[0];
rz(2.8517826) q[1];
sx q[1];
rz(-2.402014) q[1];
sx q[1];
rz(-2.0934385) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6266039) q[0];
sx q[0];
rz(-2.0467313) q[0];
sx q[0];
rz(-0.80608741) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1996908) q[2];
sx q[2];
rz(-2.2976029) q[2];
sx q[2];
rz(-2.294159) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75412616) q[1];
sx q[1];
rz(-1.4782463) q[1];
sx q[1];
rz(-2.8184163) q[1];
x q[2];
rz(-2.2506511) q[3];
sx q[3];
rz(-1.9562625) q[3];
sx q[3];
rz(2.1894574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(0.90083814) q[2];
rz(-1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(1.9074915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.95935217) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(-0.59610468) q[0];
rz(-1.6456564) q[1];
sx q[1];
rz(-0.89769617) q[1];
sx q[1];
rz(1.2449107) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0743474) q[0];
sx q[0];
rz(-1.0676358) q[0];
sx q[0];
rz(-2.7748681) q[0];
rz(0.55576022) q[2];
sx q[2];
rz(-0.55609497) q[2];
sx q[2];
rz(-0.13242002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5925496) q[1];
sx q[1];
rz(-2.2144496) q[1];
sx q[1];
rz(2.2085269) q[1];
x q[2];
rz(2.3915668) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(2.1961574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9101377) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-2.9166252) q[2];
rz(-3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(-0.47880539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(-2.8616469) q[0];
rz(1.4631368) q[1];
sx q[1];
rz(-1.8755553) q[1];
sx q[1];
rz(0.25269145) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42851105) q[0];
sx q[0];
rz(-1.6185074) q[0];
sx q[0];
rz(-1.8639355) q[0];
rz(0.018888868) q[2];
sx q[2];
rz(-0.9364555) q[2];
sx q[2];
rz(2.6332476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5340427) q[1];
sx q[1];
rz(-1.0654963) q[1];
sx q[1];
rz(2.3984548) q[1];
rz(-pi) q[2];
rz(-2.2884376) q[3];
sx q[3];
rz(-2.1907638) q[3];
sx q[3];
rz(-2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8213886) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(1.6097216) q[2];
rz(1.948471) q[3];
sx q[3];
rz(-0.90819287) q[3];
sx q[3];
rz(1.0866722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(-2.8727818) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-2.862646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19764087) q[0];
sx q[0];
rz(-0.949172) q[0];
sx q[0];
rz(0.62244121) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(-0.47889027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.61356269) q[1];
sx q[1];
rz(-2.2655305) q[1];
sx q[1];
rz(-1.0747521) q[1];
rz(-pi) q[2];
rz(-1.7989743) q[3];
sx q[3];
rz(-2.5573213) q[3];
sx q[3];
rz(0.97464558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(-1.0507978) q[3];
sx q[3];
rz(-1.2069586) q[3];
sx q[3];
rz(-1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.46362296) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(-1.8883702) q[0];
rz(0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(-1.1463096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3619986) q[0];
sx q[0];
rz(-0.75095526) q[0];
sx q[0];
rz(0.10303084) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4856604) q[2];
sx q[2];
rz(-0.99870517) q[2];
sx q[2];
rz(0.42665542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1100626) q[1];
sx q[1];
rz(-1.4974125) q[1];
sx q[1];
rz(0.88295464) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54543145) q[3];
sx q[3];
rz(-2.430393) q[3];
sx q[3];
rz(1.4382854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9877732) q[2];
sx q[2];
rz(-2.1058857) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(2.231797) q[3];
sx q[3];
rz(-2.4715021) q[3];
sx q[3];
rz(-0.36809665) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290454) q[0];
sx q[0];
rz(-0.61360252) q[0];
sx q[0];
rz(-0.01424271) q[0];
rz(-0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(-1.7600118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4603699) q[0];
sx q[0];
rz(-1.5602221) q[0];
sx q[0];
rz(-1.4282385) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6689156) q[2];
sx q[2];
rz(-2.5604904) q[2];
sx q[2];
rz(-0.89400089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8091619) q[1];
sx q[1];
rz(-2.5659487) q[1];
sx q[1];
rz(-1.3047421) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6081309) q[3];
sx q[3];
rz(-1.7603612) q[3];
sx q[3];
rz(1.61434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0649197) q[2];
sx q[2];
rz(-1.2023456) q[2];
sx q[2];
rz(0.98999611) q[2];
rz(-0.89407095) q[3];
sx q[3];
rz(-0.48864135) q[3];
sx q[3];
rz(-0.22542424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0513231) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(1.3399711) q[1];
sx q[1];
rz(-0.62146386) q[1];
sx q[1];
rz(0.38846831) q[1];
rz(1.1006533) q[2];
sx q[2];
rz(-0.91081506) q[2];
sx q[2];
rz(0.87115661) q[2];
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
