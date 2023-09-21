OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(5.8941547) q[0];
sx q[0];
rz(11.682792) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0966914) q[0];
sx q[0];
rz(-0.20151073) q[0];
sx q[0];
rz(-2.6411396) q[0];
rz(-pi) q[1];
rz(-1.2486357) q[2];
sx q[2];
rz(-0.82167168) q[2];
sx q[2];
rz(-1.8495714) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8625496) q[1];
sx q[1];
rz(-2.0837796) q[1];
sx q[1];
rz(-1.951745) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0953193) q[3];
sx q[3];
rz(-1.7414021) q[3];
sx q[3];
rz(-1.6149278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9464232) q[2];
sx q[2];
rz(-0.98313466) q[2];
sx q[2];
rz(-2.9602489) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(-0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(1.5997255) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.331859) q[0];
sx q[0];
rz(-0.74685687) q[0];
sx q[0];
rz(2.5144308) q[0];
x q[1];
rz(-1.2453169) q[2];
sx q[2];
rz(-2.6633334) q[2];
sx q[2];
rz(1.9667728) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9215556) q[1];
sx q[1];
rz(-1.6065292) q[1];
sx q[1];
rz(1.5608556) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48331355) q[3];
sx q[3];
rz(-1.1477074) q[3];
sx q[3];
rz(1.9139293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5388422) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(-1.4146457) q[2];
rz(2.291262) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(1.6833646) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26281115) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(0.61022726) q[0];
rz(-1.2894851) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-2.1496225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85342825) q[0];
sx q[0];
rz(-2.5531904) q[0];
sx q[0];
rz(-1.9073652) q[0];
x q[1];
rz(0.25408557) q[2];
sx q[2];
rz(-0.68426364) q[2];
sx q[2];
rz(2.4917045) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75533463) q[1];
sx q[1];
rz(-1.4812246) q[1];
sx q[1];
rz(-2.7002525) q[1];
rz(-pi) q[2];
rz(-1.8157585) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(-0.20416343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38301864) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(2.8733011) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85369337) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(2.7365141) q[0];
rz(-0.69008094) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-0.69782034) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989379) q[0];
sx q[0];
rz(-1.5473167) q[0];
sx q[0];
rz(1.4809181) q[0];
rz(-pi) q[1];
rz(0.58535107) q[2];
sx q[2];
rz(-0.62830892) q[2];
sx q[2];
rz(-1.5141687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88155014) q[1];
sx q[1];
rz(-1.7545732) q[1];
sx q[1];
rz(1.0764513) q[1];
rz(-pi) q[2];
rz(2.844595) q[3];
sx q[3];
rz(-1.1668491) q[3];
sx q[3];
rz(-1.7904074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.7439338) q[2];
sx q[2];
rz(1.5092124) q[2];
rz(-0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(-1.0255381) q[0];
rz(2.569596) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(-2.5122723) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48110163) q[0];
sx q[0];
rz(-1.8697303) q[0];
sx q[0];
rz(2.7116508) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25482486) q[2];
sx q[2];
rz(-2.2307768) q[2];
sx q[2];
rz(-0.13435907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3931261) q[1];
sx q[1];
rz(-1.4561635) q[1];
sx q[1];
rz(-2.6266891) q[1];
rz(-1.1961597) q[3];
sx q[3];
rz(-0.63510886) q[3];
sx q[3];
rz(2.0600968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-2.7992451) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489007) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(1.7350896) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.3669744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.470984) q[0];
sx q[0];
rz(-0.99262041) q[0];
sx q[0];
rz(-1.8355808) q[0];
rz(-pi) q[1];
rz(2.30079) q[2];
sx q[2];
rz(-1.7992939) q[2];
sx q[2];
rz(0.92323869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5806611) q[1];
sx q[1];
rz(-1.7226189) q[1];
sx q[1];
rz(0.092756943) q[1];
rz(-pi) q[2];
rz(2.572445) q[3];
sx q[3];
rz(-1.9372809) q[3];
sx q[3];
rz(1.6397427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.7657123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85526953) q[0];
sx q[0];
rz(-0.77676847) q[0];
sx q[0];
rz(-1.1458678) q[0];
x q[1];
rz(-1.6872348) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(-1.4637228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3657276) q[1];
sx q[1];
rz(-2.9850246) q[1];
sx q[1];
rz(-0.74128976) q[1];
rz(-pi) q[2];
rz(1.901058) q[3];
sx q[3];
rz(-1.5152647) q[3];
sx q[3];
rz(-0.33952573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4574796) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(-2.4462162) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59657997) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(0.73356432) q[0];
rz(0.60797524) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(-2.9072993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.308134) q[0];
sx q[0];
rz(-1.3123543) q[0];
sx q[0];
rz(-2.0636369) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2393164) q[2];
sx q[2];
rz(-0.96792816) q[2];
sx q[2];
rz(-2.5831646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56470796) q[1];
sx q[1];
rz(-1.3815834) q[1];
sx q[1];
rz(2.3343759) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1566625) q[3];
sx q[3];
rz(-1.8551991) q[3];
sx q[3];
rz(1.5581074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0155448) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(1.7162494) q[2];
rz(1.6783293) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.6459758) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5995246) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.19338) q[0];
rz(-1.880973) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(-2.1967922) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856857) q[0];
sx q[0];
rz(-0.80230306) q[0];
sx q[0];
rz(-2.5923652) q[0];
x q[1];
rz(-0.94061942) q[2];
sx q[2];
rz(-2.665328) q[2];
sx q[2];
rz(0.95048743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33465696) q[1];
sx q[1];
rz(-1.2922704) q[1];
sx q[1];
rz(1.3929699) q[1];
x q[2];
rz(-1.0778905) q[3];
sx q[3];
rz(-2.7422046) q[3];
sx q[3];
rz(2.3139017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9251359) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(-2.8533868) q[2];
rz(-2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11186803) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(2.7669725) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(2.250681) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3044395) q[0];
sx q[0];
rz(-1.3225978) q[0];
sx q[0];
rz(0.34557839) q[0];
x q[1];
rz(1.9642481) q[2];
sx q[2];
rz(-2.4911454) q[2];
sx q[2];
rz(-2.8124867) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9526457) q[1];
sx q[1];
rz(-1.9083605) q[1];
sx q[1];
rz(-1.7811716) q[1];
rz(0.37836214) q[3];
sx q[3];
rz(-2.9394657) q[3];
sx q[3];
rz(-2.9581021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23641071) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(2.6399844) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-2.0093909) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-0.51433993) q[2];
sx q[2];
rz(-0.30403501) q[2];
sx q[2];
rz(-2.4652849) q[2];
rz(3.0922079) q[3];
sx q[3];
rz(-2.1168843) q[3];
sx q[3];
rz(-2.4612853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
