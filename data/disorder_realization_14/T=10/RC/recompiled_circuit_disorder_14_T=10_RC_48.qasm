OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(0.88357893) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(4.598773) q[1];
sx q[1];
rz(7.4814319) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449013) q[0];
sx q[0];
rz(-0.20151073) q[0];
sx q[0];
rz(2.6411396) q[0];
rz(-0.3281524) q[2];
sx q[2];
rz(-0.80291623) q[2];
sx q[2];
rz(-2.3053055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4858422) q[1];
sx q[1];
rz(-1.2409004) q[1];
sx q[1];
rz(2.5962023) q[1];
x q[2];
rz(1.4000113) q[3];
sx q[3];
rz(-1.5251953) q[3];
sx q[3];
rz(-3.1053229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1951695) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-0.26120734) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3023892) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-0.38683495) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889266) q[0];
sx q[0];
rz(-0.98836556) q[0];
sx q[0];
rz(2.0684588) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16427152) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(-2.3300366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4914354) q[1];
sx q[1];
rz(-0.037089247) q[1];
sx q[1];
rz(2.8703719) q[1];
rz(-pi) q[2];
rz(-1.1003483) q[3];
sx q[3];
rz(-1.1332266) q[3];
sx q[3];
rz(-0.13089422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.26281115) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(-2.5313654) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(-2.1496225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1410071) q[0];
sx q[0];
rz(-1.7551384) q[0];
sx q[0];
rz(-2.1328451) q[0];
x q[1];
rz(0.66833468) q[2];
sx q[2];
rz(-1.7303581) q[2];
sx q[2];
rz(2.4192686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75533463) q[1];
sx q[1];
rz(-1.660368) q[1];
sx q[1];
rz(-2.7002525) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46996689) q[3];
sx q[3];
rz(-2.6372058) q[3];
sx q[3];
rz(-2.4117163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(2.8733011) q[2];
rz(-0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85369337) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-0.40507856) q[0];
rz(0.69008094) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(-0.69782034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989379) q[0];
sx q[0];
rz(-1.5942759) q[0];
sx q[0];
rz(1.6606746) q[0];
rz(-pi) q[1];
rz(-1.1890829) q[2];
sx q[2];
rz(-2.082798) q[2];
sx q[2];
rz(0.94101671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2600425) q[1];
sx q[1];
rz(-1.3870194) q[1];
sx q[1];
rz(2.0651414) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9911733) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(0.099893173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71965376) q[2];
sx q[2];
rz(-1.7439338) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-2.1160545) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(2.5122723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48110163) q[0];
sx q[0];
rz(-1.2718624) q[0];
sx q[0];
rz(-2.7116508) q[0];
rz(2.8867678) q[2];
sx q[2];
rz(-2.2307768) q[2];
sx q[2];
rz(-3.0072336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.24229901) q[1];
sx q[1];
rz(-2.0819903) q[1];
sx q[1];
rz(1.4392698) q[1];
rz(-pi) q[2];
rz(-2.878177) q[3];
sx q[3];
rz(-0.98589555) q[3];
sx q[3];
rz(1.5358621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-2.7992451) q[2];
rz(1.7957548) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(-2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(-2.956399) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-2.0506737) q[1];
sx q[1];
rz(1.7746183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6706086) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(-1.8355808) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(2.218354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5806611) q[1];
sx q[1];
rz(-1.4189737) q[1];
sx q[1];
rz(-3.0488357) q[1];
x q[2];
rz(0.6188789) q[3];
sx q[3];
rz(-2.4757909) q[3];
sx q[3];
rz(2.562079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(2.2582167) q[2];
rz(1.7287438) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(-3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.7657123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28988518) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(-0.38498199) q[0];
x q[1];
rz(2.0443633) q[2];
sx q[2];
rz(-1.5175022) q[2];
sx q[2];
rz(0.0034982861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.618305) q[1];
sx q[1];
rz(-1.6860645) q[1];
sx q[1];
rz(-1.6769888) q[1];
rz(1.7405628) q[3];
sx q[3];
rz(-0.3347291) q[3];
sx q[3];
rz(-1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6841131) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(0.69537648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59657997) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(2.5336174) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(-0.2342934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5406571) q[0];
sx q[0];
rz(-1.0957076) q[0];
sx q[0];
rz(2.8500772) q[0];
rz(2.2393164) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(-2.5831646) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5768847) q[1];
sx q[1];
rz(-1.3815834) q[1];
sx q[1];
rz(2.3343759) q[1];
rz(-2.8041744) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12604788) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(-1.4253433) q[2];
rz(-1.4632633) q[3];
sx q[3];
rz(-0.78444702) q[3];
sx q[3];
rz(1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(-1.19338) q[0];
rz(1.2606196) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017100447) q[0];
sx q[0];
rz(-1.1860577) q[0];
sx q[0];
rz(0.7229294) q[0];
rz(-pi) q[1];
rz(1.9656885) q[2];
sx q[2];
rz(-1.8443622) q[2];
sx q[2];
rz(-0.045189518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2855125) q[1];
sx q[1];
rz(-1.3998919) q[1];
sx q[1];
rz(0.28275615) q[1];
x q[2];
rz(-0.19712574) q[3];
sx q[3];
rz(-1.2211495) q[3];
sx q[3];
rz(0.2998578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21645674) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-2.8533868) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11186803) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(-2.7669725) q[1];
sx q[1];
rz(-1.7381784) q[1];
sx q[1];
rz(-0.8909117) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4961632) q[0];
sx q[0];
rz(-1.2362288) q[0];
sx q[0];
rz(1.3076925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28384039) q[2];
sx q[2];
rz(-2.1643057) q[2];
sx q[2];
rz(2.3317091) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76162321) q[1];
sx q[1];
rz(-0.39561158) q[1];
sx q[1];
rz(0.53669866) q[1];
x q[2];
rz(-0.37836214) q[3];
sx q[3];
rz(-2.9394657) q[3];
sx q[3];
rz(-0.18349056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-0.44617173) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63507737) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-0.51433993) q[2];
sx q[2];
rz(-0.30403501) q[2];
sx q[2];
rz(-2.4652849) q[2];
rz(-1.0241667) q[3];
sx q[3];
rz(-1.6129941) q[3];
sx q[3];
rz(2.2254406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
