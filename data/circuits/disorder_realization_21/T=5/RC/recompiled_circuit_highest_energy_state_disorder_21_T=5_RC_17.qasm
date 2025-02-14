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
rz(-1.9953097) q[0];
sx q[0];
rz(-0.028385552) q[0];
sx q[0];
rz(-2.1838768) q[0];
rz(-2.7081642) q[1];
sx q[1];
rz(1.537701) q[1];
sx q[1];
rz(8.796506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0032438) q[0];
sx q[0];
rz(-1.4491557) q[0];
sx q[0];
rz(-2.7697589) q[0];
x q[1];
rz(-2.6896954) q[2];
sx q[2];
rz(-2.0263645) q[2];
sx q[2];
rz(-1.6476375) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0121668) q[1];
sx q[1];
rz(-2.2291227) q[1];
sx q[1];
rz(2.7321187) q[1];
x q[2];
rz(0.099277012) q[3];
sx q[3];
rz(-1.9813207) q[3];
sx q[3];
rz(1.6308189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6750703) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(-2.7543219) q[2];
rz(-1.1747423) q[3];
sx q[3];
rz(-1.6956804) q[3];
sx q[3];
rz(-2.2430879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627477) q[0];
sx q[0];
rz(-1.400482) q[0];
sx q[0];
rz(1.8722906) q[0];
rz(-1.6954039) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(0.92099014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35996485) q[0];
sx q[0];
rz(-0.69929294) q[0];
sx q[0];
rz(2.0398606) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7537506) q[2];
sx q[2];
rz(-2.2277846) q[2];
sx q[2];
rz(-0.25568257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9780636) q[1];
sx q[1];
rz(-2.7394751) q[1];
sx q[1];
rz(-2.0744214) q[1];
rz(0.030618592) q[3];
sx q[3];
rz(-0.22600442) q[3];
sx q[3];
rz(-1.5196368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15128073) q[2];
sx q[2];
rz(-1.1996148) q[2];
sx q[2];
rz(-2.3178103) q[2];
rz(1.524205) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23564944) q[0];
sx q[0];
rz(-2.2610569) q[0];
sx q[0];
rz(-0.5916416) q[0];
rz(-0.94332424) q[1];
sx q[1];
rz(-1.0572409) q[1];
sx q[1];
rz(1.0234157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9647261) q[0];
sx q[0];
rz(-2.9826678) q[0];
sx q[0];
rz(0.83035161) q[0];
rz(-pi) q[1];
rz(2.7696174) q[2];
sx q[2];
rz(-1.0468654) q[2];
sx q[2];
rz(1.9934479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.07853807) q[1];
sx q[1];
rz(-1.5397433) q[1];
sx q[1];
rz(1.478315) q[1];
x q[2];
rz(-0.39822762) q[3];
sx q[3];
rz(-2.4169528) q[3];
sx q[3];
rz(2.6716324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25980514) q[2];
sx q[2];
rz(-1.9116348) q[2];
sx q[2];
rz(1.7477431) q[2];
rz(1.7143837) q[3];
sx q[3];
rz(-0.87427679) q[3];
sx q[3];
rz(-2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706936) q[0];
sx q[0];
rz(-0.40632668) q[0];
sx q[0];
rz(-0.64737493) q[0];
rz(-0.24230832) q[1];
sx q[1];
rz(-1.4360042) q[1];
sx q[1];
rz(1.4586331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3998916) q[0];
sx q[0];
rz(-1.8938602) q[0];
sx q[0];
rz(1.7123651) q[0];
rz(-pi) q[1];
rz(-0.83593145) q[2];
sx q[2];
rz(-2.3337337) q[2];
sx q[2];
rz(0.69617803) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9414166) q[1];
sx q[1];
rz(-1.6072909) q[1];
sx q[1];
rz(-1.0851651) q[1];
rz(-pi) q[2];
rz(2.2599873) q[3];
sx q[3];
rz(-2.2223916) q[3];
sx q[3];
rz(0.46796614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.709939) q[2];
sx q[2];
rz(-2.2971051) q[2];
sx q[2];
rz(-1.2246152) q[2];
rz(0.25105181) q[3];
sx q[3];
rz(-2.4243088) q[3];
sx q[3];
rz(-2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5421903) q[0];
sx q[0];
rz(-1.1325862) q[0];
sx q[0];
rz(-0.19790025) q[0];
rz(1.2359765) q[1];
sx q[1];
rz(-0.37207741) q[1];
sx q[1];
rz(0.0045675357) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3749927) q[0];
sx q[0];
rz(-0.99762929) q[0];
sx q[0];
rz(-1.8660924) q[0];
rz(1.4312237) q[2];
sx q[2];
rz(-2.6032631) q[2];
sx q[2];
rz(1.1280504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5089043) q[1];
sx q[1];
rz(-1.7204072) q[1];
sx q[1];
rz(-0.27638159) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7907333) q[3];
sx q[3];
rz(-1.7688732) q[3];
sx q[3];
rz(0.7367062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4602451) q[2];
sx q[2];
rz(-0.494151) q[2];
sx q[2];
rz(2.8013308) q[2];
rz(1.6081238) q[3];
sx q[3];
rz(-1.3932649) q[3];
sx q[3];
rz(1.5197915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.915864) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(1.9541784) q[0];
rz(-1.7200708) q[1];
sx q[1];
rz(-0.41821304) q[1];
sx q[1];
rz(-0.13030599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4069945) q[0];
sx q[0];
rz(-3.0049374) q[0];
sx q[0];
rz(1.1278485) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58660686) q[2];
sx q[2];
rz(-2.1663423) q[2];
sx q[2];
rz(-1.4968703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2060616) q[1];
sx q[1];
rz(-1.5701248) q[1];
sx q[1];
rz(-0.16783451) q[1];
rz(1.5339507) q[3];
sx q[3];
rz(-0.26290694) q[3];
sx q[3];
rz(-0.84883603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.80798951) q[2];
sx q[2];
rz(-1.5520059) q[2];
sx q[2];
rz(-1.0398593) q[2];
rz(2.9676159) q[3];
sx q[3];
rz(-1.1287374) q[3];
sx q[3];
rz(0.66966331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60798764) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(-0.84130353) q[0];
rz(1.3249506) q[1];
sx q[1];
rz(-2.1957896) q[1];
sx q[1];
rz(-1.673505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45638915) q[0];
sx q[0];
rz(-0.51420553) q[0];
sx q[0];
rz(-2.3715764) q[0];
x q[1];
rz(1.5505474) q[2];
sx q[2];
rz(-0.51873461) q[2];
sx q[2];
rz(-2.869019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3928396) q[1];
sx q[1];
rz(-2.3733099) q[1];
sx q[1];
rz(0.70559754) q[1];
rz(-pi) q[2];
rz(-1.2135336) q[3];
sx q[3];
rz(-1.8566086) q[3];
sx q[3];
rz(-0.88102007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2738721) q[2];
sx q[2];
rz(-0.90110675) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(0.7209512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4475187) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(-0.41282594) q[0];
rz(1.5326356) q[1];
sx q[1];
rz(-2.2282579) q[1];
sx q[1];
rz(-2.434381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21239534) q[0];
sx q[0];
rz(-1.8786504) q[0];
sx q[0];
rz(2.3677127) q[0];
rz(0.17036713) q[2];
sx q[2];
rz(-0.74935407) q[2];
sx q[2];
rz(0.10766497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0452214) q[1];
sx q[1];
rz(-1.0412843) q[1];
sx q[1];
rz(1.96223) q[1];
rz(-pi) q[2];
rz(1.8052638) q[3];
sx q[3];
rz(-0.55113652) q[3];
sx q[3];
rz(-2.7296327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1327208) q[2];
sx q[2];
rz(-0.29942313) q[2];
sx q[2];
rz(-0.51187619) q[2];
rz(0.68073186) q[3];
sx q[3];
rz(-1.3635037) q[3];
sx q[3];
rz(0.16630047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7452288) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(2.6863099) q[0];
rz(1.0633172) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(0.071970073) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0531168) q[0];
sx q[0];
rz(-1.5375111) q[0];
sx q[0];
rz(-1.5831335) q[0];
x q[1];
rz(0.44092245) q[2];
sx q[2];
rz(-1.5498232) q[2];
sx q[2];
rz(1.6787337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.94737201) q[1];
sx q[1];
rz(-2.3810219) q[1];
sx q[1];
rz(1.0494227) q[1];
rz(-pi) q[2];
rz(1.9017327) q[3];
sx q[3];
rz(-1.8171105) q[3];
sx q[3];
rz(0.34231753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7887743) q[2];
sx q[2];
rz(-2.7249551) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(-2.3173053) q[3];
sx q[3];
rz(-2.0485179) q[3];
sx q[3];
rz(2.2881959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1085827) q[0];
sx q[0];
rz(-2.0350631) q[0];
sx q[0];
rz(-0.95091096) q[0];
rz(1.7225601) q[1];
sx q[1];
rz(-2.6585237) q[1];
sx q[1];
rz(0.61666617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2717903) q[0];
sx q[0];
rz(-1.5464725) q[0];
sx q[0];
rz(-0.64016827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1918814) q[2];
sx q[2];
rz(-0.99493631) q[2];
sx q[2];
rz(2.4873231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92041278) q[1];
sx q[1];
rz(-1.3592459) q[1];
sx q[1];
rz(-1.5736363) q[1];
rz(-pi) q[2];
rz(2.1415878) q[3];
sx q[3];
rz(-2.1943008) q[3];
sx q[3];
rz(-0.68736651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0509433) q[2];
sx q[2];
rz(-1.9885149) q[2];
sx q[2];
rz(-1.0672182) q[2];
rz(-2.0459335) q[3];
sx q[3];
rz(-1.2376384) q[3];
sx q[3];
rz(-2.1777976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.511956) q[0];
sx q[0];
rz(-2.1558599) q[0];
sx q[0];
rz(2.070367) q[0];
rz(-1.4818954) q[1];
sx q[1];
rz(-2.0493458) q[1];
sx q[1];
rz(0.58072166) q[1];
rz(1.5459677) q[2];
sx q[2];
rz(-1.8332743) q[2];
sx q[2];
rz(-1.0854032) q[2];
rz(1.7076013) q[3];
sx q[3];
rz(-1.8809263) q[3];
sx q[3];
rz(3.025007) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
