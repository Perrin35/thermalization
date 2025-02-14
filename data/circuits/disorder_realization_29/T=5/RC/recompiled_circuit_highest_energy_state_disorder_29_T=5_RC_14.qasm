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
rz(-0.40773243) q[0];
sx q[0];
rz(4.3507504) q[0];
sx q[0];
rz(11.103295) q[0];
rz(0.2585803) q[1];
sx q[1];
rz(-1.9039896) q[1];
sx q[1];
rz(3.0078476) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2581078) q[0];
sx q[0];
rz(-0.15722577) q[0];
sx q[0];
rz(-1.6088465) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76579801) q[2];
sx q[2];
rz(-2.4054689) q[2];
sx q[2];
rz(0.36020261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7642149) q[1];
sx q[1];
rz(-1.0395607) q[1];
sx q[1];
rz(-1.9218535) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5695581) q[3];
sx q[3];
rz(-0.15286013) q[3];
sx q[3];
rz(0.78939702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75970834) q[2];
sx q[2];
rz(-2.3298161) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(-1.367502) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1644208) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(-1.0154065) q[0];
rz(2.7482391) q[1];
sx q[1];
rz(-1.4721847) q[1];
sx q[1];
rz(-0.70158395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4751101) q[0];
sx q[0];
rz(-2.283264) q[0];
sx q[0];
rz(3.0300167) q[0];
rz(1.267822) q[2];
sx q[2];
rz(-1.5605443) q[2];
sx q[2];
rz(1.645719) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1893411) q[1];
sx q[1];
rz(-0.32508141) q[1];
sx q[1];
rz(-1.3334683) q[1];
rz(-2.3823678) q[3];
sx q[3];
rz(-2.1442778) q[3];
sx q[3];
rz(-2.218117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.027448805) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(1.1055498) q[2];
rz(0.28087273) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(2.9893379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0356782) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(-2.0689082) q[0];
rz(3.1283123) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(-2.5650011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7130797) q[0];
sx q[0];
rz(-1.858485) q[0];
sx q[0];
rz(-1.6384533) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0835952) q[2];
sx q[2];
rz(-1.4382677) q[2];
sx q[2];
rz(-0.26187632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62200786) q[1];
sx q[1];
rz(-1.7652006) q[1];
sx q[1];
rz(-1.0854875) q[1];
x q[2];
rz(-1.9535096) q[3];
sx q[3];
rz(-0.59186223) q[3];
sx q[3];
rz(-2.446021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2836634) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(0.64209765) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-0.27500209) q[3];
sx q[3];
rz(-0.16998418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2791486) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(-0.24236648) q[0];
rz(0.53388059) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(1.4885611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72172644) q[0];
sx q[0];
rz(-0.97630608) q[0];
sx q[0];
rz(-1.6564903) q[0];
rz(1.2514417) q[2];
sx q[2];
rz(-1.9803626) q[2];
sx q[2];
rz(-1.3495654) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3900534) q[1];
sx q[1];
rz(-1.281847) q[1];
sx q[1];
rz(2.3759936) q[1];
x q[2];
rz(-2.8418301) q[3];
sx q[3];
rz(-1.6090172) q[3];
sx q[3];
rz(2.770407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36294595) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(1.6710336) q[2];
rz(3.0342024) q[3];
sx q[3];
rz(-1.3067747) q[3];
sx q[3];
rz(1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.030468) q[0];
sx q[0];
rz(-0.46115369) q[0];
sx q[0];
rz(1.2029458) q[0];
rz(2.2466834) q[1];
sx q[1];
rz(-1.9590961) q[1];
sx q[1];
rz(0.21200171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1269147) q[0];
sx q[0];
rz(-2.0058332) q[0];
sx q[0];
rz(2.2382909) q[0];
x q[1];
rz(-2.6190125) q[2];
sx q[2];
rz(-1.8327288) q[2];
sx q[2];
rz(-0.22522989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3983237) q[1];
sx q[1];
rz(-1.124427) q[1];
sx q[1];
rz(2.2012987) q[1];
rz(1.0174973) q[3];
sx q[3];
rz(-1.5626255) q[3];
sx q[3];
rz(3.0752575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8733069) q[2];
sx q[2];
rz(-0.77750677) q[2];
sx q[2];
rz(-0.11534616) q[2];
rz(-0.98571944) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91144052) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(-2.4440515) q[0];
rz(0.41360924) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(-2.4381309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0516104) q[0];
sx q[0];
rz(-1.4282266) q[0];
sx q[0];
rz(0.27326126) q[0];
x q[1];
rz(1.037961) q[2];
sx q[2];
rz(-0.67801266) q[2];
sx q[2];
rz(-0.17282669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78402987) q[1];
sx q[1];
rz(-2.5270487) q[1];
sx q[1];
rz(-2.4667506) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88075068) q[3];
sx q[3];
rz(-0.38057571) q[3];
sx q[3];
rz(0.56169034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0387663) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(1.9290257) q[2];
rz(-2.8596527) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(-0.50659242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7956227) q[0];
sx q[0];
rz(-2.3891734) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(2.1961191) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(2.4623154) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69791171) q[0];
sx q[0];
rz(-2.6583932) q[0];
sx q[0];
rz(-2.9761821) q[0];
x q[1];
rz(-1.8273583) q[2];
sx q[2];
rz(-0.96254331) q[2];
sx q[2];
rz(1.7781374) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0019466) q[1];
sx q[1];
rz(-1.0602385) q[1];
sx q[1];
rz(-2.9278048) q[1];
rz(-pi) q[2];
rz(0.019370989) q[3];
sx q[3];
rz(-1.8550411) q[3];
sx q[3];
rz(-1.5694048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2519553) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(1.4944705) q[2];
rz(2.5877118) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.7249829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81903356) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(-1.092528) q[0];
rz(0.93387261) q[1];
sx q[1];
rz(-1.4051508) q[1];
sx q[1];
rz(-0.78528231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0781495) q[0];
sx q[0];
rz(-0.64029303) q[0];
sx q[0];
rz(-2.3904353) q[0];
rz(2.2319517) q[2];
sx q[2];
rz(-2.4850436) q[2];
sx q[2];
rz(0.54443371) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31639659) q[1];
sx q[1];
rz(-1.8421208) q[1];
sx q[1];
rz(-1.3153988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1340895) q[3];
sx q[3];
rz(-2.1434382) q[3];
sx q[3];
rz(0.2084032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2698722) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(2.5093057) q[2];
rz(-0.91442433) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(0.44347611) q[0];
rz(0.30287287) q[1];
sx q[1];
rz(-0.58285204) q[1];
sx q[1];
rz(1.833896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46638495) q[0];
sx q[0];
rz(-1.8716011) q[0];
sx q[0];
rz(1.7487026) q[0];
x q[1];
rz(1.2595749) q[2];
sx q[2];
rz(-1.0838795) q[2];
sx q[2];
rz(-2.4608909) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5622632) q[1];
sx q[1];
rz(-2.0585357) q[1];
sx q[1];
rz(-0.30178782) q[1];
x q[2];
rz(-2.1122718) q[3];
sx q[3];
rz(-0.53526894) q[3];
sx q[3];
rz(-2.1524371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0206454) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(-0.81765085) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(-0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870975) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(-2.7823271) q[0];
rz(0.66419762) q[1];
sx q[1];
rz(-1.8134873) q[1];
sx q[1];
rz(2.7204303) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065124113) q[0];
sx q[0];
rz(-1.9389922) q[0];
sx q[0];
rz(2.9028331) q[0];
x q[1];
rz(-3.1362651) q[2];
sx q[2];
rz(-1.6432646) q[2];
sx q[2];
rz(0.42124149) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1637472) q[1];
sx q[1];
rz(-2.6085269) q[1];
sx q[1];
rz(-2.1043489) q[1];
x q[2];
rz(0.75712689) q[3];
sx q[3];
rz(-2.3806678) q[3];
sx q[3];
rz(-1.8426551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8757214) q[2];
sx q[2];
rz(-0.096780626) q[2];
sx q[2];
rz(-1.2145915) q[2];
rz(0.38861361) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(1.0987561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6471967) q[0];
sx q[0];
rz(-1.5955441) q[0];
sx q[0];
rz(1.2936976) q[0];
rz(2.9248059) q[1];
sx q[1];
rz(-0.68991359) q[1];
sx q[1];
rz(0.51295113) q[1];
rz(-2.0016128) q[2];
sx q[2];
rz(-0.93968334) q[2];
sx q[2];
rz(0.46776007) q[2];
rz(0.46617266) q[3];
sx q[3];
rz(-1.8099347) q[3];
sx q[3];
rz(3.1113887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
