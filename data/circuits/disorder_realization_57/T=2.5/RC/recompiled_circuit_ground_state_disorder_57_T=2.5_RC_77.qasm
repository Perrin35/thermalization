OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2630513) q[0];
sx q[0];
rz(3.75293) q[0];
sx q[0];
rz(11.073025) q[0];
rz(-2.2892294) q[1];
sx q[1];
rz(-1.747793) q[1];
sx q[1];
rz(-1.9115023) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6540007) q[0];
sx q[0];
rz(-2.0474259) q[0];
sx q[0];
rz(-2.0311902) q[0];
rz(0.70582055) q[2];
sx q[2];
rz(-2.5748944) q[2];
sx q[2];
rz(1.3659878) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7609357) q[1];
sx q[1];
rz(-2.0208997) q[1];
sx q[1];
rz(-0.37827079) q[1];
x q[2];
rz(1.2870077) q[3];
sx q[3];
rz(-1.5012245) q[3];
sx q[3];
rz(-0.33659914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7386231) q[2];
sx q[2];
rz(-1.7454742) q[2];
sx q[2];
rz(1.1633066) q[2];
rz(-2.2033384) q[3];
sx q[3];
rz(-0.80621755) q[3];
sx q[3];
rz(1.4279385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52629483) q[0];
sx q[0];
rz(-1.2328923) q[0];
sx q[0];
rz(1.1616608) q[0];
rz(1.7192526) q[1];
sx q[1];
rz(-1.4931581) q[1];
sx q[1];
rz(3.0551522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0491701) q[0];
sx q[0];
rz(-1.0455275) q[0];
sx q[0];
rz(-0.18245936) q[0];
rz(-pi) q[1];
rz(-2.6506054) q[2];
sx q[2];
rz(-1.7684968) q[2];
sx q[2];
rz(-0.035015496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.045043) q[1];
sx q[1];
rz(-0.031686671) q[1];
sx q[1];
rz(-2.6331054) q[1];
x q[2];
rz(-1.9374886) q[3];
sx q[3];
rz(-0.18309284) q[3];
sx q[3];
rz(-2.2568373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57261673) q[2];
sx q[2];
rz(-2.6628351) q[2];
sx q[2];
rz(-3.0968481) q[2];
rz(-2.4217126) q[3];
sx q[3];
rz(-1.5174815) q[3];
sx q[3];
rz(-1.450479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2429263) q[0];
sx q[0];
rz(-1.1761605) q[0];
sx q[0];
rz(1.9157008) q[0];
rz(2.2236845) q[1];
sx q[1];
rz(-1.2914912) q[1];
sx q[1];
rz(-1.8570522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96439904) q[0];
sx q[0];
rz(-1.5209235) q[0];
sx q[0];
rz(2.0727488) q[0];
rz(-pi) q[1];
x q[1];
rz(1.578737) q[2];
sx q[2];
rz(-2.1924332) q[2];
sx q[2];
rz(1.4162743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0513902) q[1];
sx q[1];
rz(-1.0388684) q[1];
sx q[1];
rz(0.73511413) q[1];
rz(-pi) q[2];
rz(-0.35090943) q[3];
sx q[3];
rz(-2.4882462) q[3];
sx q[3];
rz(-0.86337435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43117493) q[2];
sx q[2];
rz(-1.147889) q[2];
sx q[2];
rz(-2.5244024) q[2];
rz(0.91444531) q[3];
sx q[3];
rz(-0.72943288) q[3];
sx q[3];
rz(2.7845553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451839) q[0];
sx q[0];
rz(-0.04627385) q[0];
sx q[0];
rz(-2.1531877) q[0];
rz(1.4811966) q[1];
sx q[1];
rz(-2.587187) q[1];
sx q[1];
rz(0.60715094) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6364215) q[0];
sx q[0];
rz(-1.8431516) q[0];
sx q[0];
rz(0.26557458) q[0];
rz(-pi) q[1];
rz(-0.49272196) q[2];
sx q[2];
rz(-0.70198503) q[2];
sx q[2];
rz(-1.5240325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.41736351) q[1];
sx q[1];
rz(-2.0290142) q[1];
sx q[1];
rz(0.37727892) q[1];
rz(3.0444381) q[3];
sx q[3];
rz(-1.3911671) q[3];
sx q[3];
rz(0.068795242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0561169) q[2];
sx q[2];
rz(-1.3885219) q[2];
sx q[2];
rz(0.2307387) q[2];
rz(-0.76656109) q[3];
sx q[3];
rz(-2.9975588) q[3];
sx q[3];
rz(1.8421596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6799927) q[0];
sx q[0];
rz(-1.6277286) q[0];
sx q[0];
rz(-0.067721279) q[0];
rz(1.8374775) q[1];
sx q[1];
rz(-1.291052) q[1];
sx q[1];
rz(-1.4363323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3910113) q[0];
sx q[0];
rz(-0.89871797) q[0];
sx q[0];
rz(-0.4990905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7737232) q[2];
sx q[2];
rz(-2.6393916) q[2];
sx q[2];
rz(1.7562438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7244514) q[1];
sx q[1];
rz(-1.0994689) q[1];
sx q[1];
rz(0.36496867) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9201305) q[3];
sx q[3];
rz(-1.2173183) q[3];
sx q[3];
rz(2.2686951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0453673) q[2];
sx q[2];
rz(-2.2787091) q[2];
sx q[2];
rz(0.53675845) q[2];
rz(1.9219575) q[3];
sx q[3];
rz(-0.9044956) q[3];
sx q[3];
rz(1.1902683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14559513) q[0];
sx q[0];
rz(-2.1639316) q[0];
sx q[0];
rz(-1.4033432) q[0];
rz(0.504269) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(2.9020342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1288236) q[0];
sx q[0];
rz(-1.4137303) q[0];
sx q[0];
rz(1.7536909) q[0];
rz(-pi) q[1];
rz(2.7272322) q[2];
sx q[2];
rz(-1.0584497) q[2];
sx q[2];
rz(3.1238585) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7799222) q[1];
sx q[1];
rz(-2.4677271) q[1];
sx q[1];
rz(3.0008631) q[1];
x q[2];
rz(0.40294874) q[3];
sx q[3];
rz(-1.1004173) q[3];
sx q[3];
rz(1.7960447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2704894) q[2];
sx q[2];
rz(-1.5364545) q[2];
sx q[2];
rz(2.8883873) q[2];
rz(2.1829055) q[3];
sx q[3];
rz(-2.2264693) q[3];
sx q[3];
rz(1.2747214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7764928) q[0];
sx q[0];
rz(-2.4204142) q[0];
sx q[0];
rz(1.391063) q[0];
rz(-0.59025383) q[1];
sx q[1];
rz(-1.9275459) q[1];
sx q[1];
rz(0.50277695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8322635) q[0];
sx q[0];
rz(-1.676275) q[0];
sx q[0];
rz(3.1270887) q[0];
rz(-1.9421135) q[2];
sx q[2];
rz(-0.94778143) q[2];
sx q[2];
rz(-2.3065904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8233899) q[1];
sx q[1];
rz(-2.3651337) q[1];
sx q[1];
rz(1.9253325) q[1];
x q[2];
rz(1.8842949) q[3];
sx q[3];
rz(-1.0491706) q[3];
sx q[3];
rz(-0.23594638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34660029) q[2];
sx q[2];
rz(-1.9575926) q[2];
sx q[2];
rz(0.45197519) q[2];
rz(0.31451264) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(0.048129169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5101584) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(1.6360224) q[0];
rz(-0.10558852) q[1];
sx q[1];
rz(-2.0629203) q[1];
sx q[1];
rz(-0.67965913) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9980695) q[0];
sx q[0];
rz(-0.98833246) q[0];
sx q[0];
rz(0.67654718) q[0];
x q[1];
rz(2.7397268) q[2];
sx q[2];
rz(-0.99152126) q[2];
sx q[2];
rz(1.4039672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7235106) q[1];
sx q[1];
rz(-2.6716304) q[1];
sx q[1];
rz(-1.577792) q[1];
rz(-pi) q[2];
rz(0.18631492) q[3];
sx q[3];
rz(-1.3491231) q[3];
sx q[3];
rz(1.9566389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6203561) q[2];
sx q[2];
rz(-2.607589) q[2];
sx q[2];
rz(-0.83068577) q[2];
rz(2.1031117) q[3];
sx q[3];
rz(-0.65282789) q[3];
sx q[3];
rz(1.4000819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164383) q[0];
sx q[0];
rz(-1.9898299) q[0];
sx q[0];
rz(1.0333767) q[0];
rz(2.0694464) q[1];
sx q[1];
rz(-1.9629581) q[1];
sx q[1];
rz(2.5796366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009223) q[0];
sx q[0];
rz(-1.7096247) q[0];
sx q[0];
rz(0.85900659) q[0];
rz(-pi) q[1];
rz(0.91479723) q[2];
sx q[2];
rz(-1.2848228) q[2];
sx q[2];
rz(-0.15483072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6889557) q[1];
sx q[1];
rz(-1.9105043) q[1];
sx q[1];
rz(2.3477786) q[1];
rz(-3.1329145) q[3];
sx q[3];
rz(-1.5132705) q[3];
sx q[3];
rz(0.044139095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.098103913) q[2];
sx q[2];
rz(-1.6489112) q[2];
sx q[2];
rz(1.7334422) q[2];
rz(-0.060338542) q[3];
sx q[3];
rz(-1.8294168) q[3];
sx q[3];
rz(-1.8797125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9982346) q[0];
sx q[0];
rz(-2.3682605) q[0];
sx q[0];
rz(-1.7711357) q[0];
rz(-1.0073608) q[1];
sx q[1];
rz(-1.7528563) q[1];
sx q[1];
rz(0.14152424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62011224) q[0];
sx q[0];
rz(-2.4864288) q[0];
sx q[0];
rz(-0.76458365) q[0];
rz(-pi) q[1];
rz(-2.3958489) q[2];
sx q[2];
rz(-0.088938449) q[2];
sx q[2];
rz(-2.3740785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87301577) q[1];
sx q[1];
rz(-2.0076224) q[1];
sx q[1];
rz(-0.42307968) q[1];
rz(2.0662289) q[3];
sx q[3];
rz(-1.5513707) q[3];
sx q[3];
rz(2.7904449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0274028) q[2];
sx q[2];
rz(-1.9582615) q[2];
sx q[2];
rz(1.2664504) q[2];
rz(3.0053084) q[3];
sx q[3];
rz(-0.91629052) q[3];
sx q[3];
rz(-0.19792476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65205735) q[0];
sx q[0];
rz(-1.0132402) q[0];
sx q[0];
rz(-1.5860438) q[0];
rz(2.1736705) q[1];
sx q[1];
rz(-0.5236917) q[1];
sx q[1];
rz(2.0655469) q[1];
rz(-0.81071557) q[2];
sx q[2];
rz(-1.492112) q[2];
sx q[2];
rz(0.44853609) q[2];
rz(-0.051911671) q[3];
sx q[3];
rz(-1.7604699) q[3];
sx q[3];
rz(-1.5552675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
