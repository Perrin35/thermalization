OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(-2.5728777) q[0];
rz(1.2110127) q[1];
sx q[1];
rz(-2.14415) q[1];
sx q[1];
rz(2.8740191) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5934323) q[0];
sx q[0];
rz(-2.555937) q[0];
sx q[0];
rz(-0.36958739) q[0];
rz(-1.9047584) q[2];
sx q[2];
rz(-1.6777473) q[2];
sx q[2];
rz(3.1283875) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5423746) q[1];
sx q[1];
rz(-1.5732048) q[1];
sx q[1];
rz(1.5791248) q[1];
x q[2];
rz(1.7115643) q[3];
sx q[3];
rz(-1.632431) q[3];
sx q[3];
rz(0.34256645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.25207511) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(1.1313103) q[2];
rz(0.45025292) q[3];
sx q[3];
rz(-0.69142747) q[3];
sx q[3];
rz(2.1972307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.4985519) q[0];
sx q[0];
rz(-0.93198553) q[0];
sx q[0];
rz(-1.7339647) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.718037) q[1];
sx q[1];
rz(2.5462467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4099497) q[0];
sx q[0];
rz(-0.34706693) q[0];
sx q[0];
rz(2.3763229) q[0];
rz(-pi) q[1];
rz(-0.96244241) q[2];
sx q[2];
rz(-2.3909937) q[2];
sx q[2];
rz(-0.19718328) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7253984) q[1];
sx q[1];
rz(-0.42896118) q[1];
sx q[1];
rz(3.129175) q[1];
x q[2];
rz(1.1717623) q[3];
sx q[3];
rz(-2.1431987) q[3];
sx q[3];
rz(0.31203285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88196102) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(2.372443) q[2];
rz(1.7800356) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(-0.59534016) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24460569) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(-2.8047674) q[0];
rz(1.7103851) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(-3.0444042) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6568156) q[0];
sx q[0];
rz(-2.6069399) q[0];
sx q[0];
rz(1.9770245) q[0];
x q[1];
rz(0.23938208) q[2];
sx q[2];
rz(-0.39153831) q[2];
sx q[2];
rz(-1.7042314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1657451) q[1];
sx q[1];
rz(-0.37384181) q[1];
sx q[1];
rz(-0.38623078) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52690701) q[3];
sx q[3];
rz(-0.95427536) q[3];
sx q[3];
rz(-1.34672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97239697) q[2];
sx q[2];
rz(-1.7652067) q[2];
sx q[2];
rz(0.81673679) q[2];
rz(-0.9225325) q[3];
sx q[3];
rz(-2.7042992) q[3];
sx q[3];
rz(1.1923265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836477) q[0];
sx q[0];
rz(-1.8035996) q[0];
sx q[0];
rz(-0.86679593) q[0];
rz(-1.2358933) q[1];
sx q[1];
rz(-2.0265323) q[1];
sx q[1];
rz(-0.24196504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45296882) q[0];
sx q[0];
rz(-1.3039661) q[0];
sx q[0];
rz(2.2950606) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1461805) q[2];
sx q[2];
rz(-0.17689366) q[2];
sx q[2];
rz(1.5695656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6852808) q[1];
sx q[1];
rz(-1.9748828) q[1];
sx q[1];
rz(-1.7194952) q[1];
rz(-pi) q[2];
rz(-0.091598467) q[3];
sx q[3];
rz(-2.6428416) q[3];
sx q[3];
rz(2.2745511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26291004) q[2];
sx q[2];
rz(-1.6308558) q[2];
sx q[2];
rz(-2.6061457) q[2];
rz(0.49324909) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(0.44805995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0754452) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(2.2976663) q[0];
rz(1.6663724) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(2.7484238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12411453) q[0];
sx q[0];
rz(-2.3459166) q[0];
sx q[0];
rz(1.9301729) q[0];
rz(-pi) q[1];
rz(-1.7885782) q[2];
sx q[2];
rz(-2.1960207) q[2];
sx q[2];
rz(-0.54131258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0991126) q[1];
sx q[1];
rz(-2.5951324) q[1];
sx q[1];
rz(2.5853755) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.353998) q[3];
sx q[3];
rz(-1.9301842) q[3];
sx q[3];
rz(-1.5415292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7273442) q[2];
sx q[2];
rz(-0.20432893) q[2];
sx q[2];
rz(-2.7434529) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-0.78854338) q[3];
sx q[3];
rz(-1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2274461) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(2.2858802) q[0];
rz(-0.69264597) q[1];
sx q[1];
rz(-2.1470862) q[1];
sx q[1];
rz(-2.8725502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97996431) q[0];
sx q[0];
rz(-1.3792896) q[0];
sx q[0];
rz(-3.1167517) q[0];
rz(-1.6330209) q[2];
sx q[2];
rz(-2.6341558) q[2];
sx q[2];
rz(2.1560046) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4743691) q[1];
sx q[1];
rz(-1.3880236) q[1];
sx q[1];
rz(-2.3007042) q[1];
x q[2];
rz(-2.2773197) q[3];
sx q[3];
rz(-0.24989299) q[3];
sx q[3];
rz(2.8003729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0272224) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(3.031292) q[2];
rz(-0.86841622) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(-1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39349839) q[0];
sx q[0];
rz(-2.6415934) q[0];
sx q[0];
rz(-0.86135832) q[0];
rz(-1.512108) q[1];
sx q[1];
rz(-1.4605099) q[1];
sx q[1];
rz(0.82383627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.829946) q[0];
sx q[0];
rz(-0.59893805) q[0];
sx q[0];
rz(-1.7700559) q[0];
x q[1];
rz(-1.1608949) q[2];
sx q[2];
rz(-0.86386743) q[2];
sx q[2];
rz(0.17690578) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2715335) q[1];
sx q[1];
rz(-0.64836577) q[1];
sx q[1];
rz(-0.35761498) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7384375) q[3];
sx q[3];
rz(-2.4488827) q[3];
sx q[3];
rz(1.5628536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6350101) q[2];
sx q[2];
rz(-2.625605) q[2];
sx q[2];
rz(-2.3487976) q[2];
rz(3.1363764) q[3];
sx q[3];
rz(-2.3514533) q[3];
sx q[3];
rz(-1.4504356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.97013) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(-2.9520853) q[0];
rz(-2.3686523) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(2.5453087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7211821) q[0];
sx q[0];
rz(-0.89119688) q[0];
sx q[0];
rz(1.3534989) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79769602) q[2];
sx q[2];
rz(-1.0587278) q[2];
sx q[2];
rz(1.860581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3459732) q[1];
sx q[1];
rz(-1.6633777) q[1];
sx q[1];
rz(3.0365192) q[1];
x q[2];
rz(0.52395405) q[3];
sx q[3];
rz(-0.92902459) q[3];
sx q[3];
rz(0.061836035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4160055) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(1.2131946) q[2];
rz(-0.98617918) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(-1.2590316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1213433) q[0];
sx q[0];
rz(-1.5556524) q[0];
sx q[0];
rz(1.1100618) q[0];
rz(0.32866651) q[1];
sx q[1];
rz(-1.5866491) q[1];
sx q[1];
rz(1.8448255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27889869) q[0];
sx q[0];
rz(-1.4771013) q[0];
sx q[0];
rz(0.91178943) q[0];
rz(-pi) q[1];
rz(1.718156) q[2];
sx q[2];
rz(-0.40032101) q[2];
sx q[2];
rz(-1.4555228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83767855) q[1];
sx q[1];
rz(-2.0635567) q[1];
sx q[1];
rz(-1.4999092) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7862919) q[3];
sx q[3];
rz(-2.2249376) q[3];
sx q[3];
rz(0.69132016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0729596) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(0.10406058) q[2];
rz(2.0067298) q[3];
sx q[3];
rz(-1.7770504) q[3];
sx q[3];
rz(2.9141736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65570152) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(-2.4110598) q[0];
rz(0.55745521) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(-2.7117859) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107501) q[0];
sx q[0];
rz(-1.5702899) q[0];
sx q[0];
rz(0.7617801) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5556106) q[2];
sx q[2];
rz(-0.44011099) q[2];
sx q[2];
rz(-0.42428478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.046172) q[1];
sx q[1];
rz(-0.37994994) q[1];
sx q[1];
rz(-1.3543966) q[1];
x q[2];
rz(-1.881322) q[3];
sx q[3];
rz(-0.38703296) q[3];
sx q[3];
rz(-2.455445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3320015) q[2];
sx q[2];
rz(-1.0217228) q[2];
sx q[2];
rz(2.8487955) q[2];
rz(3.0012567) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(-0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705227) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(-0.71939214) q[1];
sx q[1];
rz(-1.8665301) q[1];
sx q[1];
rz(-2.9449609) q[1];
rz(-1.0631845) q[2];
sx q[2];
rz(-1.7287935) q[2];
sx q[2];
rz(-2.4274735) q[2];
rz(3.0075913) q[3];
sx q[3];
rz(-1.8905427) q[3];
sx q[3];
rz(-0.34656634) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
