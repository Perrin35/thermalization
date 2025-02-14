OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(3.1867653) q[0];
sx q[0];
rz(10.09633) q[0];
rz(2.1454732) q[1];
sx q[1];
rz(5.4944333) q[1];
sx q[1];
rz(6.5715437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.45123) q[0];
sx q[0];
rz(-1.7685862) q[0];
sx q[0];
rz(3.1044699) q[0];
rz(-pi) q[1];
rz(0.34121193) q[2];
sx q[2];
rz(-0.69319572) q[2];
sx q[2];
rz(-2.5643519) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25569281) q[1];
sx q[1];
rz(-1.5330557) q[1];
sx q[1];
rz(-0.44817145) q[1];
rz(1.7156473) q[3];
sx q[3];
rz(-1.1205691) q[3];
sx q[3];
rz(-2.1046771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4396189) q[2];
sx q[2];
rz(-2.7957323) q[2];
sx q[2];
rz(2.4227179) q[2];
rz(-1.4465205) q[3];
sx q[3];
rz(-1.4868163) q[3];
sx q[3];
rz(-2.1046624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0816536) q[0];
sx q[0];
rz(-2.7103598) q[0];
sx q[0];
rz(-2.8531139) q[0];
rz(0.55229315) q[1];
sx q[1];
rz(-2.0497132) q[1];
sx q[1];
rz(-1.9658032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73574725) q[0];
sx q[0];
rz(-1.7205392) q[0];
sx q[0];
rz(-1.2558054) q[0];
x q[1];
rz(-0.98227882) q[2];
sx q[2];
rz(-0.20766307) q[2];
sx q[2];
rz(1.3192624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1604483) q[1];
sx q[1];
rz(-1.4466043) q[1];
sx q[1];
rz(0.9915413) q[1];
rz(1.5351686) q[3];
sx q[3];
rz(-0.032101121) q[3];
sx q[3];
rz(0.42313448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4484078) q[2];
sx q[2];
rz(-1.2625445) q[2];
sx q[2];
rz(3.1191678) q[2];
rz(-1.4261931) q[3];
sx q[3];
rz(-1.8636999) q[3];
sx q[3];
rz(0.9001596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4876323) q[0];
sx q[0];
rz(-1.910169) q[0];
sx q[0];
rz(2.3768429) q[0];
rz(1.359831) q[1];
sx q[1];
rz(-1.2218852) q[1];
sx q[1];
rz(-1.8870032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0645197) q[0];
sx q[0];
rz(-1.4409587) q[0];
sx q[0];
rz(-3.0373851) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7255035) q[2];
sx q[2];
rz(-0.35230428) q[2];
sx q[2];
rz(-1.2809629) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4855328) q[1];
sx q[1];
rz(-1.1120218) q[1];
sx q[1];
rz(-1.8937673) q[1];
x q[2];
rz(2.7990667) q[3];
sx q[3];
rz(-1.4493128) q[3];
sx q[3];
rz(-3.1241724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4059056) q[2];
sx q[2];
rz(-0.76852208) q[2];
sx q[2];
rz(-3.1206257) q[2];
rz(1.5603125) q[3];
sx q[3];
rz(-1.0617278) q[3];
sx q[3];
rz(0.90644065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2727994) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(-1.772076) q[0];
rz(2.4319793) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(-2.4443464) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0727834) q[0];
sx q[0];
rz(-1.8801196) q[0];
sx q[0];
rz(0.42865045) q[0];
x q[1];
rz(-2.4575649) q[2];
sx q[2];
rz(-2.1993756) q[2];
sx q[2];
rz(0.39234871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62909758) q[1];
sx q[1];
rz(-1.8726908) q[1];
sx q[1];
rz(-0.90934335) q[1];
x q[2];
rz(1.3342821) q[3];
sx q[3];
rz(-1.4035038) q[3];
sx q[3];
rz(3.0574556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-2.1944025) q[2];
sx q[2];
rz(0.67898018) q[2];
rz(-2.0164356) q[3];
sx q[3];
rz(-1.9966639) q[3];
sx q[3];
rz(-2.9185435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515902) q[0];
sx q[0];
rz(-1.9380049) q[0];
sx q[0];
rz(-0.077127174) q[0];
rz(1.9902825) q[1];
sx q[1];
rz(-1.5733893) q[1];
sx q[1];
rz(0.77879771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4753195) q[0];
sx q[0];
rz(-0.059564807) q[0];
sx q[0];
rz(0.93713974) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16806099) q[2];
sx q[2];
rz(-1.2894783) q[2];
sx q[2];
rz(2.9945053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8647789) q[1];
sx q[1];
rz(-1.9052231) q[1];
sx q[1];
rz(-0.2590551) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058249931) q[3];
sx q[3];
rz(-0.16213972) q[3];
sx q[3];
rz(-1.8193965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6246346) q[2];
sx q[2];
rz(-1.7004852) q[2];
sx q[2];
rz(-1.1901633) q[2];
rz(0.55365753) q[3];
sx q[3];
rz(-0.62479574) q[3];
sx q[3];
rz(-0.73738086) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5364285) q[0];
sx q[0];
rz(-1.1506511) q[0];
sx q[0];
rz(-0.27106699) q[0];
rz(-1.0003264) q[1];
sx q[1];
rz(-1.9258291) q[1];
sx q[1];
rz(0.86885524) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2216744) q[0];
sx q[0];
rz(-2.1561047) q[0];
sx q[0];
rz(-1.0006389) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0539114) q[2];
sx q[2];
rz(-0.73747915) q[2];
sx q[2];
rz(-0.35754851) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7160401) q[1];
sx q[1];
rz(-1.2231266) q[1];
sx q[1];
rz(2.9613609) q[1];
rz(-pi) q[2];
rz(1.1017076) q[3];
sx q[3];
rz(-0.84381754) q[3];
sx q[3];
rz(-0.25694822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.038593682) q[2];
sx q[2];
rz(-0.30240348) q[2];
sx q[2];
rz(-1.137286) q[2];
rz(-2.8310827) q[3];
sx q[3];
rz(-0.91028428) q[3];
sx q[3];
rz(-2.212132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1167269) q[0];
sx q[0];
rz(-1.7127345) q[0];
sx q[0];
rz(-1.0799991) q[0];
rz(2.179821) q[1];
sx q[1];
rz(-2.5764143) q[1];
sx q[1];
rz(-2.9186509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7540383) q[0];
sx q[0];
rz(-1.6056653) q[0];
sx q[0];
rz(3.1341482) q[0];
rz(-pi) q[1];
rz(-0.30855631) q[2];
sx q[2];
rz(-0.94538222) q[2];
sx q[2];
rz(0.61185123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60513105) q[1];
sx q[1];
rz(-2.2830182) q[1];
sx q[1];
rz(0.3682767) q[1];
x q[2];
rz(0.97496521) q[3];
sx q[3];
rz(-1.5830073) q[3];
sx q[3];
rz(1.7705767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1792533) q[2];
sx q[2];
rz(-0.27955678) q[2];
sx q[2];
rz(1.9635828) q[2];
rz(-2.8152605) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(-1.1403722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2278263) q[0];
sx q[0];
rz(-0.95781177) q[0];
sx q[0];
rz(-2.3334099) q[0];
rz(-2.783964) q[1];
sx q[1];
rz(-1.459815) q[1];
sx q[1];
rz(2.0590032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21511997) q[0];
sx q[0];
rz(-1.400504) q[0];
sx q[0];
rz(-2.2567458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9679734) q[2];
sx q[2];
rz(-0.589314) q[2];
sx q[2];
rz(0.1620503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33629234) q[1];
sx q[1];
rz(-2.3545579) q[1];
sx q[1];
rz(1.8563104) q[1];
x q[2];
rz(2.0478422) q[3];
sx q[3];
rz(-1.1594605) q[3];
sx q[3];
rz(0.75440948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.064934405) q[2];
sx q[2];
rz(-0.4883464) q[2];
sx q[2];
rz(1.2154382) q[2];
rz(0.91228929) q[3];
sx q[3];
rz(-0.89022294) q[3];
sx q[3];
rz(-0.54273763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8488345) q[0];
sx q[0];
rz(-2.500535) q[0];
sx q[0];
rz(0.42375281) q[0];
rz(-1.9641701) q[1];
sx q[1];
rz(-2.2260428) q[1];
sx q[1];
rz(-2.7659069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9388449) q[0];
sx q[0];
rz(-1.3934982) q[0];
sx q[0];
rz(3.0696654) q[0];
x q[1];
rz(0.10924364) q[2];
sx q[2];
rz(-2.901361) q[2];
sx q[2];
rz(1.1775398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2089952) q[1];
sx q[1];
rz(-1.0699341) q[1];
sx q[1];
rz(-0.074525699) q[1];
x q[2];
rz(1.9104092) q[3];
sx q[3];
rz(-2.2582128) q[3];
sx q[3];
rz(2.825045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11289135) q[2];
sx q[2];
rz(-1.4914923) q[2];
sx q[2];
rz(-2.4493307) q[2];
rz(2.4568457) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(-0.14028604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5570062) q[0];
sx q[0];
rz(-0.96730119) q[0];
sx q[0];
rz(1.6242356) q[0];
rz(-1.5380305) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(1.921152) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.206911) q[0];
sx q[0];
rz(-0.62536541) q[0];
sx q[0];
rz(-1.4689133) q[0];
rz(-pi) q[1];
rz(-0.50165711) q[2];
sx q[2];
rz(-1.3634472) q[2];
sx q[2];
rz(-0.32688552) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0076589) q[1];
sx q[1];
rz(-0.87532212) q[1];
sx q[1];
rz(1.201406) q[1];
rz(-2.4046005) q[3];
sx q[3];
rz(-1.3496282) q[3];
sx q[3];
rz(-1.2545409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5662235) q[2];
sx q[2];
rz(-1.0074002) q[2];
sx q[2];
rz(2.2929906) q[2];
rz(0.97950116) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(-0.037467329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10850567) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(-0.55950821) q[1];
sx q[1];
rz(-0.63897501) q[1];
sx q[1];
rz(-0.41313304) q[1];
rz(1.826936) q[2];
sx q[2];
rz(-1.8825681) q[2];
sx q[2];
rz(-2.067461) q[2];
rz(-0.091690334) q[3];
sx q[3];
rz(-2.5602362) q[3];
sx q[3];
rz(0.14701281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
