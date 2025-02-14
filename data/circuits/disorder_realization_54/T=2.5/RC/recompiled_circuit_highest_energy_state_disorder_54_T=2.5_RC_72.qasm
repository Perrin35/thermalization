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
rz(2.8132791) q[0];
sx q[0];
rz(-2.0067196) q[0];
sx q[0];
rz(2.5757134) q[0];
rz(0.24428754) q[1];
sx q[1];
rz(-1.3574418) q[1];
sx q[1];
rz(0.53425962) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7843648) q[0];
sx q[0];
rz(-2.2450232) q[0];
sx q[0];
rz(2.7927464) q[0];
x q[1];
rz(2.7083394) q[2];
sx q[2];
rz(-0.64068551) q[2];
sx q[2];
rz(-2.6034466) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45785832) q[1];
sx q[1];
rz(-2.8401655) q[1];
sx q[1];
rz(-3.0617759) q[1];
rz(-pi) q[2];
rz(-2.684428) q[3];
sx q[3];
rz(-1.7948397) q[3];
sx q[3];
rz(0.30686298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1285706) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(-2.8924083) q[2];
rz(1.3439882) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(-0.070076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.2263612) q[0];
sx q[0];
rz(2.1517854) q[0];
rz(2.3509332) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(-0.31165037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383331) q[0];
sx q[0];
rz(-2.0796392) q[0];
sx q[0];
rz(2.3340387) q[0];
rz(1.3840394) q[2];
sx q[2];
rz(-1.4644564) q[2];
sx q[2];
rz(-2.3723536) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6271882) q[1];
sx q[1];
rz(-1.7429105) q[1];
sx q[1];
rz(1.8132416) q[1];
x q[2];
rz(0.95590638) q[3];
sx q[3];
rz(-1.8528189) q[3];
sx q[3];
rz(2.7789555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.1602465) q[2];
sx q[2];
rz(2.1902093) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-0.97672668) q[3];
sx q[3];
rz(1.2736646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-3.0271295) q[0];
rz(-2.139367) q[1];
sx q[1];
rz(-2.7148235) q[1];
sx q[1];
rz(0.1618298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.072542) q[0];
sx q[0];
rz(-0.74063535) q[0];
sx q[0];
rz(-2.1669751) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6587018) q[2];
sx q[2];
rz(-1.1832596) q[2];
sx q[2];
rz(1.0532925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56282957) q[1];
sx q[1];
rz(-0.96841267) q[1];
sx q[1];
rz(-1.8627326) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1597685) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(-2.4730554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74730045) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(3.0008345) q[2];
rz(0.55177871) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9363339) q[0];
sx q[0];
rz(-0.2660428) q[0];
sx q[0];
rz(2.4667013) q[0];
rz(-2.9871125) q[1];
sx q[1];
rz(-1.4090425) q[1];
sx q[1];
rz(2.6836269) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7688528) q[0];
sx q[0];
rz(-2.495666) q[0];
sx q[0];
rz(2.4619308) q[0];
x q[1];
rz(0.30720023) q[2];
sx q[2];
rz(-1.9824902) q[2];
sx q[2];
rz(-1.9780985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6667216) q[1];
sx q[1];
rz(-2.7242766) q[1];
sx q[1];
rz(2.6732207) q[1];
rz(-pi) q[2];
rz(-2.9513957) q[3];
sx q[3];
rz(-0.51589291) q[3];
sx q[3];
rz(0.35515337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0080879) q[2];
sx q[2];
rz(-0.68024457) q[2];
sx q[2];
rz(2.8221455) q[2];
rz(-1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1595681) q[0];
sx q[0];
rz(-2.0052795) q[0];
sx q[0];
rz(-1.9200448) q[0];
rz(2.9087861) q[1];
sx q[1];
rz(-2.5722952) q[1];
sx q[1];
rz(0.98308841) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1637835) q[0];
sx q[0];
rz(-1.0484694) q[0];
sx q[0];
rz(-2.1420826) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4120502) q[2];
sx q[2];
rz(-0.6266784) q[2];
sx q[2];
rz(1.145663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5447588) q[1];
sx q[1];
rz(-1.7167517) q[1];
sx q[1];
rz(2.9533141) q[1];
x q[2];
rz(1.2125254) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(3.0606396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7588707) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(2.7663084) q[2];
rz(1.075607) q[3];
sx q[3];
rz(-2.7552102) q[3];
sx q[3];
rz(2.6066499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4009092) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(-1.4373454) q[0];
rz(0.078723343) q[1];
sx q[1];
rz(-1.7648141) q[1];
sx q[1];
rz(0.54814235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054436) q[0];
sx q[0];
rz(-0.6983122) q[0];
sx q[0];
rz(0.33238073) q[0];
rz(0.55193837) q[2];
sx q[2];
rz(-1.9935529) q[2];
sx q[2];
rz(1.1814337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74433148) q[1];
sx q[1];
rz(-1.9588699) q[1];
sx q[1];
rz(-2.4281888) q[1];
x q[2];
rz(-2.1918212) q[3];
sx q[3];
rz(-2.2249537) q[3];
sx q[3];
rz(-2.014924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4569725) q[2];
sx q[2];
rz(-2.520884) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(1.5843377) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-0.61711446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.5439344) q[0];
sx q[0];
rz(-1.1827396) q[0];
sx q[0];
rz(0.4959929) q[0];
rz(-1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(1.1167663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3163213) q[0];
sx q[0];
rz(-2.6351235) q[0];
sx q[0];
rz(-0.78790085) q[0];
rz(-1.1363813) q[2];
sx q[2];
rz(-0.58131733) q[2];
sx q[2];
rz(-2.8535247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3254814) q[1];
sx q[1];
rz(-1.8891836) q[1];
sx q[1];
rz(0.35811425) q[1];
rz(0.12567606) q[3];
sx q[3];
rz(-0.45914859) q[3];
sx q[3];
rz(2.7530376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.83069673) q[2];
sx q[2];
rz(-2.1572957) q[2];
sx q[2];
rz(0.35487077) q[2];
rz(-2.3762459) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4527721) q[0];
sx q[0];
rz(-1.6625762) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(2.2701524) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(-1.8188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3324748) q[0];
sx q[0];
rz(-1.5189384) q[0];
sx q[0];
rz(2.8807667) q[0];
rz(-pi) q[1];
rz(2.9051612) q[2];
sx q[2];
rz(-2.2364103) q[2];
sx q[2];
rz(-0.98275358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8812528) q[1];
sx q[1];
rz(-1.5543206) q[1];
sx q[1];
rz(2.6332098) q[1];
rz(-pi) q[2];
rz(1.7962097) q[3];
sx q[3];
rz(-1.7724049) q[3];
sx q[3];
rz(-0.20568337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4837997) q[2];
sx q[2];
rz(-1.0427661) q[2];
sx q[2];
rz(-2.9198809) q[2];
rz(2.0486369) q[3];
sx q[3];
rz(-1.7287798) q[3];
sx q[3];
rz(-2.4634585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(0.27716032) q[0];
rz(-2.3932338) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(1.1433196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4410981) q[0];
sx q[0];
rz(-1.3702533) q[0];
sx q[0];
rz(-1.2593377) q[0];
rz(-pi) q[1];
rz(0.21094811) q[2];
sx q[2];
rz(-0.8767414) q[2];
sx q[2];
rz(-1.2488493) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4797693) q[1];
sx q[1];
rz(-2.7854438) q[1];
sx q[1];
rz(2.0355862) q[1];
rz(-pi) q[2];
rz(-1.3995029) q[3];
sx q[3];
rz(-1.2138324) q[3];
sx q[3];
rz(-1.5689284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2271759) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(3.0926256) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-2.5132096) q[3];
sx q[3];
rz(-0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12298909) q[0];
sx q[0];
rz(-1.6139231) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(0.77876577) q[1];
sx q[1];
rz(-0.6929144) q[1];
sx q[1];
rz(-1.5787517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7514141) q[0];
sx q[0];
rz(-1.7008243) q[0];
sx q[0];
rz(3.0517654) q[0];
rz(0.18264211) q[2];
sx q[2];
rz(-0.96871766) q[2];
sx q[2];
rz(-1.4878359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73247468) q[1];
sx q[1];
rz(-2.1555087) q[1];
sx q[1];
rz(0.52701108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7354129) q[3];
sx q[3];
rz(-0.74062982) q[3];
sx q[3];
rz(0.1333789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(0.16555244) q[2];
rz(2.5340269) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3792569) q[0];
sx q[0];
rz(-0.48185928) q[0];
sx q[0];
rz(3.0961105) q[0];
rz(-1.634585) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(0.14590895) q[2];
sx q[2];
rz(-2.3608531) q[2];
sx q[2];
rz(0.94028801) q[2];
rz(0.60579469) q[3];
sx q[3];
rz(-1.7346564) q[3];
sx q[3];
rz(0.32614542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
