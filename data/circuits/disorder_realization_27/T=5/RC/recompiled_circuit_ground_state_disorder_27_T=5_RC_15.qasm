OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7967427) q[0];
sx q[0];
rz(-2.8673708) q[0];
sx q[0];
rz(2.5728777) q[0];
rz(-1.93058) q[1];
sx q[1];
rz(-0.99744263) q[1];
sx q[1];
rz(-2.8740191) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9833004) q[0];
sx q[0];
rz(-2.112297) q[0];
sx q[0];
rz(1.805961) q[0];
rz(-pi) q[1];
rz(1.2368343) q[2];
sx q[2];
rz(-1.6777473) q[2];
sx q[2];
rz(3.1283875) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59921801) q[1];
sx q[1];
rz(-1.5683878) q[1];
sx q[1];
rz(-1.5624678) q[1];
x q[2];
rz(-1.7115643) q[3];
sx q[3];
rz(-1.632431) q[3];
sx q[3];
rz(-0.34256645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25207511) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(2.0102823) q[2];
rz(-0.45025292) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(-0.94436193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(1.7339647) q[0];
rz(0.44250008) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(0.59534591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93564088) q[0];
sx q[0];
rz(-1.3229587) q[0];
sx q[0];
rz(-1.3252844) q[0];
x q[1];
rz(0.91752618) q[2];
sx q[2];
rz(-1.1703614) q[2];
sx q[2];
rz(-1.296907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9756979) q[1];
sx q[1];
rz(-1.5656316) q[1];
sx q[1];
rz(-0.42893202) q[1];
rz(-pi) q[2];
rz(-1.9698304) q[3];
sx q[3];
rz(-0.99839393) q[3];
sx q[3];
rz(2.8295598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.88196102) q[2];
sx q[2];
rz(-2.5233848) q[2];
sx q[2];
rz(0.76914966) q[2];
rz(1.361557) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(0.59534016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460569) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(2.8047674) q[0];
rz(-1.7103851) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(-0.097188458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48477706) q[0];
sx q[0];
rz(-2.6069399) q[0];
sx q[0];
rz(1.9770245) q[0];
rz(2.7601542) q[2];
sx q[2];
rz(-1.4801916) q[2];
sx q[2];
rz(0.355313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5775902) q[1];
sx q[1];
rz(-1.9158984) q[1];
sx q[1];
rz(-1.4240828) q[1];
rz(0.95366565) q[3];
sx q[3];
rz(-0.78804555) q[3];
sx q[3];
rz(-2.5833481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97239697) q[2];
sx q[2];
rz(-1.376386) q[2];
sx q[2];
rz(0.81673679) q[2];
rz(-2.2190602) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(1.1923265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(2.8996276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075027) q[0];
sx q[0];
rz(-0.76341141) q[0];
sx q[0];
rz(-1.1795189) q[0];
rz(-pi) q[1];
rz(1.9954122) q[2];
sx q[2];
rz(-2.964699) q[2];
sx q[2];
rz(1.572027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0556524) q[1];
sx q[1];
rz(-1.4341518) q[1];
sx q[1];
rz(-0.40811347) q[1];
x q[2];
rz(1.5210152) q[3];
sx q[3];
rz(-2.0672653) q[3];
sx q[3];
rz(-2.1703326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8786826) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(-2.6061457) q[2];
rz(2.6483436) q[3];
sx q[3];
rz(-0.89087629) q[3];
sx q[3];
rz(-2.6935327) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0754452) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(0.84392631) q[0];
rz(-1.6663724) q[1];
sx q[1];
rz(-1.6553144) q[1];
sx q[1];
rz(2.7484238) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12411453) q[0];
sx q[0];
rz(-0.79567608) q[0];
sx q[0];
rz(-1.2114197) q[0];
x q[1];
rz(-2.5049823) q[2];
sx q[2];
rz(-1.746897) q[2];
sx q[2];
rz(-2.2409093) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.04248) q[1];
sx q[1];
rz(-2.5951324) q[1];
sx q[1];
rz(-2.5853755) q[1];
rz(2.6541071) q[3];
sx q[3];
rz(-0.84934399) q[3];
sx q[3];
rz(2.8340428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7273442) q[2];
sx q[2];
rz(-0.20432893) q[2];
sx q[2];
rz(-0.3981398) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-0.78854338) q[3];
sx q[3];
rz(-1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(-0.8557125) q[0];
rz(-2.4489467) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(-2.8725502) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1616283) q[0];
sx q[0];
rz(-1.3792896) q[0];
sx q[0];
rz(-3.1167517) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6330209) q[2];
sx q[2];
rz(-2.6341558) q[2];
sx q[2];
rz(0.98558805) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4384296) q[1];
sx q[1];
rz(-0.7483349) q[1];
sx q[1];
rz(-1.8412043) q[1];
rz(-pi) q[2];
rz(1.7625436) q[3];
sx q[3];
rz(-1.4095528) q[3];
sx q[3];
rz(1.2210326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1143703) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(-3.031292) q[2];
rz(0.86841622) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(-1.4364012) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39349839) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(0.86135832) q[0];
rz(-1.512108) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(2.3177564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.829946) q[0];
sx q[0];
rz(-0.59893805) q[0];
sx q[0];
rz(1.3715368) q[0];
x q[1];
rz(-1.9806978) q[2];
sx q[2];
rz(-2.2777252) q[2];
sx q[2];
rz(0.17690578) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7303402) q[1];
sx q[1];
rz(-1.7837875) q[1];
sx q[1];
rz(2.5243882) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4895913) q[3];
sx q[3];
rz(-1.8240415) q[3];
sx q[3];
rz(2.8324003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6350101) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(-2.3487976) q[2];
rz(-3.1363764) q[3];
sx q[3];
rz(-2.3514533) q[3];
sx q[3];
rz(1.4504356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1714627) q[0];
sx q[0];
rz(-1.1928394) q[0];
sx q[0];
rz(-2.9520853) q[0];
rz(0.7729404) q[1];
sx q[1];
rz(-0.49630061) q[1];
sx q[1];
rz(-2.5453087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8533289) q[0];
sx q[0];
rz(-1.402308) q[0];
sx q[0];
rz(-0.69126076) q[0];
x q[1];
rz(2.3438966) q[2];
sx q[2];
rz(-1.0587278) q[2];
sx q[2];
rz(1.860581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6365123) q[1];
sx q[1];
rz(-0.13992913) q[1];
sx q[1];
rz(-0.72461463) q[1];
x q[2];
rz(0.98084456) q[3];
sx q[3];
rz(-2.337237) q[3];
sx q[3];
rz(0.70589069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72558713) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(1.928398) q[2];
rz(0.98617918) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(-1.8825611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0202494) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(2.0315309) q[0];
rz(0.32866651) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(1.2967671) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9220306) q[0];
sx q[0];
rz(-2.2264105) q[0];
sx q[0];
rz(3.0232885) q[0];
rz(-1.9672212) q[2];
sx q[2];
rz(-1.5135445) q[2];
sx q[2];
rz(-3.1210085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83767855) q[1];
sx q[1];
rz(-1.078036) q[1];
sx q[1];
rz(-1.6416835) q[1];
rz(-pi) q[2];
rz(1.7862919) q[3];
sx q[3];
rz(-2.2249376) q[3];
sx q[3];
rz(-0.69132016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0729596) q[2];
sx q[2];
rz(-0.98560846) q[2];
sx q[2];
rz(-3.0375321) q[2];
rz(1.1348628) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(-0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65570152) q[0];
sx q[0];
rz(-2.1567397) q[0];
sx q[0];
rz(-0.73053288) q[0];
rz(2.5841374) q[1];
sx q[1];
rz(-1.971222) q[1];
sx q[1];
rz(-2.7117859) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9011079) q[0];
sx q[0];
rz(-0.76178023) q[0];
sx q[0];
rz(-0.00073379993) q[0];
rz(-pi) q[1];
rz(-2.5556106) q[2];
sx q[2];
rz(-0.44011099) q[2];
sx q[2];
rz(0.42428478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.046172) q[1];
sx q[1];
rz(-0.37994994) q[1];
sx q[1];
rz(-1.3543966) q[1];
rz(-pi) q[2];
rz(1.2602706) q[3];
sx q[3];
rz(-2.7545597) q[3];
sx q[3];
rz(-0.68614764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80959117) q[2];
sx q[2];
rz(-1.0217228) q[2];
sx q[2];
rz(-2.8487955) q[2];
rz(0.14033595) q[3];
sx q[3];
rz(-1.0293181) q[3];
sx q[3];
rz(-0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.77107) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(2.4222005) q[1];
sx q[1];
rz(-1.8665301) q[1];
sx q[1];
rz(-2.9449609) q[1];
rz(-2.9612598) q[2];
sx q[2];
rz(-1.0700995) q[2];
sx q[2];
rz(2.3722103) q[2];
rz(0.1340014) q[3];
sx q[3];
rz(-1.25105) q[3];
sx q[3];
rz(2.7950263) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
