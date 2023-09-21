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
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6057974) q[0];
sx q[0];
rz(-1.3942766) q[0];
sx q[0];
rz(-1.4730886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3281524) q[2];
sx q[2];
rz(-2.3386764) q[2];
sx q[2];
rz(-2.3053055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6557505) q[1];
sx q[1];
rz(-1.2409004) q[1];
sx q[1];
rz(0.54539036) q[1];
rz(-pi) q[2];
rz(1.3084859) q[3];
sx q[3];
rz(-0.17671083) q[3];
sx q[3];
rz(-1.7929329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(2.9602489) q[2];
rz(-0.26120734) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-2.3852824) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5418672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5889266) q[0];
sx q[0];
rz(-0.98836556) q[0];
sx q[0];
rz(1.0731339) q[0];
rz(-1.1142271) q[2];
sx q[2];
rz(-1.7184988) q[2];
sx q[2];
rz(-2.4545124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6501573) q[1];
sx q[1];
rz(-3.1045034) q[1];
sx q[1];
rz(0.27122072) q[1];
rz(-2.0412444) q[3];
sx q[3];
rz(-1.1332266) q[3];
sx q[3];
rz(0.13089422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60275045) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(-1.7269469) q[2];
rz(2.291262) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787815) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(-1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-0.99197018) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2881644) q[0];
sx q[0];
rz(-2.5531904) q[0];
sx q[0];
rz(1.2342274) q[0];
x q[1];
rz(-1.3685554) q[2];
sx q[2];
rz(-0.91245203) q[2];
sx q[2];
rz(0.973268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62832824) q[1];
sx q[1];
rz(-0.4497512) q[1];
sx q[1];
rz(-2.9343534) q[1];
rz(-1.8157585) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(-0.20416343) q[3];
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
rz(-0.26829159) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(-2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85369337) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(2.4515117) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(0.69782034) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989379) q[0];
sx q[0];
rz(-1.5473167) q[0];
sx q[0];
rz(-1.6606746) q[0];
rz(-1.1890829) q[2];
sx q[2];
rz(-1.0587947) q[2];
sx q[2];
rz(-0.94101671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2600425) q[1];
sx q[1];
rz(-1.7545732) q[1];
sx q[1];
rz(-2.0651414) q[1];
rz(-pi) q[2];
rz(-2.1711369) q[3];
sx q[3];
rz(-2.6451023) q[3];
sx q[3];
rz(-1.1288201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(2.737282) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25800911) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-2.1160545) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(-0.62932032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2239265) q[0];
sx q[0];
rz(-1.1610982) q[0];
sx q[0];
rz(-1.8976589) q[0];
rz(0.25482486) q[2];
sx q[2];
rz(-0.91081589) q[2];
sx q[2];
rz(-3.0072336) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.021796062) q[1];
sx q[1];
rz(-2.6152059) q[1];
sx q[1];
rz(-2.9119133) q[1];
rz(-2.1719645) q[3];
sx q[3];
rz(-1.789635) q[3];
sx q[3];
rz(0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(2.7992451) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(0.19601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.3669744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.470984) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(1.3060119) q[0];
rz(0.84080266) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(-0.92323869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0095014) q[1];
sx q[1];
rz(-2.9638634) q[1];
sx q[1];
rz(1.0264261) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6188789) q[3];
sx q[3];
rz(-2.4757909) q[3];
sx q[3];
rz(-2.562079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(-2.2582167) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25154034) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(-1.863377) q[0];
rz(0.037840769) q[1];
sx q[1];
rz(-0.81532878) q[1];
sx q[1];
rz(-1.7657123) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28988518) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(0.38498199) q[0];
rz(-3.0817229) q[2];
sx q[2];
rz(-2.0436358) q[2];
sx q[2];
rz(1.5470031) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5232877) q[1];
sx q[1];
rz(-1.6860645) q[1];
sx q[1];
rz(-1.6769888) q[1];
rz(-pi) q[2];
rz(1.4010299) q[3];
sx q[3];
rz(-0.3347291) q[3];
sx q[3];
rz(1.7498121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(-0.73105556) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59657997) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(2.5336174) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(-2.9072993) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5406571) q[0];
sx q[0];
rz(-1.0957076) q[0];
sx q[0];
rz(0.29151543) q[0];
rz(0.72006165) q[2];
sx q[2];
rz(-1.0349627) q[2];
sx q[2];
rz(-0.59130397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5768847) q[1];
sx q[1];
rz(-1.7600093) q[1];
sx q[1];
rz(2.3343759) q[1];
rz(1.0844564) q[3];
sx q[3];
rz(-2.497695) q[3];
sx q[3];
rz(0.38734303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.12604788) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(1.4253433) q[2];
rz(-1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5995246) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(-1.19338) q[0];
rz(-1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9076618) q[0];
sx q[0];
rz(-2.230927) q[0];
sx q[0];
rz(2.0659167) q[0];
rz(0.94061942) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(0.95048743) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2855125) q[1];
sx q[1];
rz(-1.7417007) q[1];
sx q[1];
rz(2.8588365) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0637022) q[3];
sx q[3];
rz(-0.3993881) q[3];
sx q[3];
rz(-2.3139017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(2.8533868) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297246) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(-0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(-0.8909117) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64542949) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(-1.3076925) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1773445) q[2];
sx q[2];
rz(-2.4911454) q[2];
sx q[2];
rz(2.8124867) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3799694) q[1];
sx q[1];
rz(-0.39561158) q[1];
sx q[1];
rz(-2.604894) q[1];
rz(-pi) q[2];
rz(1.4952412) q[3];
sx q[3];
rz(-1.3831426) q[3];
sx q[3];
rz(2.5726266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(-0.50160828) q[2];
rz(1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-1.1322017) q[1];
sx q[1];
rz(-2.3846346) q[1];
sx q[1];
rz(0.089288575) q[1];
rz(1.7239465) q[2];
sx q[2];
rz(-1.3071123) q[2];
sx q[2];
rz(1.2109962) q[2];
rz(1.0241667) q[3];
sx q[3];
rz(-1.5285986) q[3];
sx q[3];
rz(-0.91615208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
