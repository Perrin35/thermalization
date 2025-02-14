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
rz(2.7338602) q[0];
sx q[0];
rz(-1.2091577) q[0];
sx q[0];
rz(-1.6785167) q[0];
rz(3.4001729) q[1];
sx q[1];
rz(1.9039896) q[1];
sx q[1];
rz(9.2910329) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2966327) q[0];
sx q[0];
rz(-1.7279074) q[0];
sx q[0];
rz(0.0060307302) q[0];
rz(-pi) q[1];
rz(2.5630579) q[2];
sx q[2];
rz(-1.0867439) q[2];
sx q[2];
rz(2.5501188) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37686302) q[1];
sx q[1];
rz(-1.2697744) q[1];
sx q[1];
rz(2.5824598) q[1];
x q[2];
rz(-0.00019076983) q[3];
sx q[3];
rz(-1.7236563) q[3];
sx q[3];
rz(2.3534485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75970834) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(-1.7740907) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(-0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1644208) q[0];
sx q[0];
rz(-1.9855969) q[0];
sx q[0];
rz(2.1261862) q[0];
rz(-2.7482391) q[1];
sx q[1];
rz(-1.4721847) q[1];
sx q[1];
rz(0.70158395) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.836228) q[0];
sx q[0];
rz(-0.71963632) q[0];
sx q[0];
rz(-1.4426065) q[0];
x q[1];
rz(-1.267822) q[2];
sx q[2];
rz(-1.5605443) q[2];
sx q[2];
rz(1.4958737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2021705) q[1];
sx q[1];
rz(-1.2551474) q[1];
sx q[1];
rz(-3.062518) q[1];
rz(-2.3823678) q[3];
sx q[3];
rz(-0.99731481) q[3];
sx q[3];
rz(-0.92347565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1141438) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(1.1055498) q[2];
rz(2.8607199) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(-2.9893379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1059145) q[0];
sx q[0];
rz(-0.67258251) q[0];
sx q[0];
rz(-2.0689082) q[0];
rz(-0.013280344) q[1];
sx q[1];
rz(-1.5417136) q[1];
sx q[1];
rz(2.5650011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.018533) q[0];
sx q[0];
rz(-1.5059239) q[0];
sx q[0];
rz(-2.8532802) q[0];
rz(-pi) q[1];
rz(-3.0579975) q[2];
sx q[2];
rz(-1.4382677) q[2];
sx q[2];
rz(-2.8797163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84725897) q[1];
sx q[1];
rz(-1.0953961) q[1];
sx q[1];
rz(-0.21902276) q[1];
rz(-pi) q[2];
rz(-1.9535096) q[3];
sx q[3];
rz(-2.5497304) q[3];
sx q[3];
rz(2.446021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85792929) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(-0.64209765) q[2];
rz(0.56504956) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-0.16998418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86244407) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(2.8992262) q[0];
rz(-0.53388059) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(-1.4885611) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3405995) q[0];
sx q[0];
rz(-1.4998319) q[0];
sx q[0];
rz(2.5453955) q[0];
x q[1];
rz(-0.62612957) q[2];
sx q[2];
rz(-2.6278547) q[2];
sx q[2];
rz(-1.0990248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0341314) q[1];
sx q[1];
rz(-2.3337769) q[1];
sx q[1];
rz(-0.40523578) q[1];
rz(-pi) q[2];
rz(-2.8418301) q[3];
sx q[3];
rz(-1.6090172) q[3];
sx q[3];
rz(2.770407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36294595) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(1.6710336) q[2];
rz(-0.10739022) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(1.86514) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11112467) q[0];
sx q[0];
rz(-0.46115369) q[0];
sx q[0];
rz(1.2029458) q[0];
rz(-0.89490926) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(-0.21200171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1231736) q[0];
sx q[0];
rz(-2.1668262) q[0];
sx q[0];
rz(-2.6072748) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49290979) q[2];
sx q[2];
rz(-0.57905902) q[2];
sx q[2];
rz(0.92307239) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0088259) q[1];
sx q[1];
rz(-1.0100875) q[1];
sx q[1];
rz(-2.6067023) q[1];
x q[2];
rz(-0.0096036951) q[3];
sx q[3];
rz(-2.1240747) q[3];
sx q[3];
rz(-1.6320848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8733069) q[2];
sx q[2];
rz(-0.77750677) q[2];
sx q[2];
rz(0.11534616) q[2];
rz(-0.98571944) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(-0.43578291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91144052) q[0];
sx q[0];
rz(-1.4608915) q[0];
sx q[0];
rz(-0.69754115) q[0];
rz(-0.41360924) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(2.4381309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1915775) q[0];
sx q[0];
rz(-2.8342025) q[0];
sx q[0];
rz(-0.48883518) q[0];
rz(-pi) q[1];
rz(-2.1773019) q[2];
sx q[2];
rz(-1.8950771) q[2];
sx q[2];
rz(-1.8285268) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3575628) q[1];
sx q[1];
rz(-0.61454391) q[1];
sx q[1];
rz(-2.4667506) q[1];
rz(-2.892214) q[3];
sx q[3];
rz(-1.2802534) q[3];
sx q[3];
rz(0.16502608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10282639) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(-1.9290257) q[2];
rz(2.8596527) q[3];
sx q[3];
rz(-1.5423256) q[3];
sx q[3];
rz(2.6350002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3459699) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(0.94547358) q[1];
sx q[1];
rz(-1.5739601) q[1];
sx q[1];
rz(2.4623154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0196456) q[0];
sx q[0];
rz(-1.4942193) q[0];
sx q[0];
rz(-0.47756739) q[0];
rz(1.3142344) q[2];
sx q[2];
rz(-2.1790493) q[2];
sx q[2];
rz(-1.7781374) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1396461) q[1];
sx q[1];
rz(-2.0813542) q[1];
sx q[1];
rz(2.9278048) q[1];
x q[2];
rz(-1.5045937) q[3];
sx q[3];
rz(-0.28488628) q[3];
sx q[3];
rz(-1.5004304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2519553) q[2];
sx q[2];
rz(-0.95459443) q[2];
sx q[2];
rz(1.6471222) q[2];
rz(-0.55388081) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.7249829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3225591) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(-2.0490647) q[0];
rz(0.93387261) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(-2.3563103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781495) q[0];
sx q[0];
rz(-0.64029303) q[0];
sx q[0];
rz(-2.3904353) q[0];
rz(-pi) q[1];
rz(-0.90964093) q[2];
sx q[2];
rz(-2.4850436) q[2];
sx q[2];
rz(-2.5971589) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.685647) q[1];
sx q[1];
rz(-0.37044493) q[1];
sx q[1];
rz(2.4042995) q[1];
rz(-2.5232072) q[3];
sx q[3];
rz(-1.2073642) q[3];
sx q[3];
rz(1.6101162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87172047) q[2];
sx q[2];
rz(-1.6951268) q[2];
sx q[2];
rz(-2.5093057) q[2];
rz(-2.2271683) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(2.9818592) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81894994) q[0];
sx q[0];
rz(-2.2845848) q[0];
sx q[0];
rz(2.6981165) q[0];
rz(0.30287287) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(1.3076967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752077) q[0];
sx q[0];
rz(-1.8716011) q[0];
sx q[0];
rz(-1.3928901) q[0];
x q[1];
rz(1.2595749) q[2];
sx q[2];
rz(-2.0577132) q[2];
sx q[2];
rz(2.4608909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.16567) q[1];
sx q[1];
rz(-2.5745086) q[1];
sx q[1];
rz(-1.0600863) q[1];
rz(-pi) q[2];
rz(2.0409706) q[3];
sx q[3];
rz(-1.8368145) q[3];
sx q[3];
rz(-0.10423743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12094721) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(2.3239418) q[3];
sx q[3];
rz(-2.652467) q[3];
sx q[3];
rz(-2.8221655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15449512) q[0];
sx q[0];
rz(-1.3767865) q[0];
sx q[0];
rz(2.7823271) q[0];
rz(-2.477395) q[1];
sx q[1];
rz(-1.8134873) q[1];
sx q[1];
rz(-0.42116234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065124113) q[0];
sx q[0];
rz(-1.2026005) q[0];
sx q[0];
rz(2.9028331) q[0];
x q[1];
rz(1.6440516) q[2];
sx q[2];
rz(-0.072663531) q[2];
sx q[2];
rz(0.49468985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56257427) q[1];
sx q[1];
rz(-2.0236602) q[1];
sx q[1];
rz(0.29154204) q[1];
rz(-2.1499885) q[3];
sx q[3];
rz(-2.0957855) q[3];
sx q[3];
rz(-0.92574173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26587129) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(-1.2145915) q[2];
rz(-2.752979) q[3];
sx q[3];
rz(-2.2190084) q[3];
sx q[3];
rz(2.0428366) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.494396) q[0];
sx q[0];
rz(-1.5955441) q[0];
sx q[0];
rz(1.2936976) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(-0.67736652) q[2];
sx q[2];
rz(-1.9147041) q[2];
sx q[2];
rz(-0.83821581) q[2];
rz(-1.8372336) q[3];
sx q[3];
rz(-2.022701) q[3];
sx q[3];
rz(-1.4823784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
