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
rz(1.1462829) q[0];
sx q[0];
rz(-3.1132071) q[0];
sx q[0];
rz(-0.95771587) q[0];
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
rz(1.4074898) q[0];
sx q[0];
rz(-0.39034319) q[0];
sx q[0];
rz(2.8170308) q[0];
rz(-pi) q[1];
rz(-2.6896954) q[2];
sx q[2];
rz(-2.0263645) q[2];
sx q[2];
rz(1.4939552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5124197) q[1];
sx q[1];
rz(-0.75885526) q[1];
sx q[1];
rz(-1.0954129) q[1];
rz(-1.3468993) q[3];
sx q[3];
rz(-0.42169398) q[3];
sx q[3];
rz(1.8754011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6750703) q[2];
sx q[2];
rz(-1.0202946) q[2];
sx q[2];
rz(-0.38727078) q[2];
rz(1.9668503) q[3];
sx q[3];
rz(-1.6956804) q[3];
sx q[3];
rz(-2.2430879) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3627477) q[0];
sx q[0];
rz(-1.400482) q[0];
sx q[0];
rz(1.8722906) q[0];
rz(1.4461888) q[1];
sx q[1];
rz(-1.8480999) q[1];
sx q[1];
rz(0.92099014) q[1];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.3878421) q[2];
sx q[2];
rz(-2.2277846) q[2];
sx q[2];
rz(0.25568257) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9780636) q[1];
sx q[1];
rz(-0.40211758) q[1];
sx q[1];
rz(-2.0744214) q[1];
rz(-pi) q[2];
rz(2.9156906) q[3];
sx q[3];
rz(-1.5639362) q[3];
sx q[3];
rz(-0.021319162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15128073) q[2];
sx q[2];
rz(-1.1996148) q[2];
sx q[2];
rz(0.82378236) q[2];
rz(1.6173877) q[3];
sx q[3];
rz(-1.5533605) q[3];
sx q[3];
rz(-1.2113781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23564944) q[0];
sx q[0];
rz(-2.2610569) q[0];
sx q[0];
rz(2.549951) q[0];
rz(-0.94332424) q[1];
sx q[1];
rz(-2.0843518) q[1];
sx q[1];
rz(-1.0234157) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34020326) q[0];
sx q[0];
rz(-1.6777628) q[0];
sx q[0];
rz(-1.4530327) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1259427) q[2];
sx q[2];
rz(-1.890939) q[2];
sx q[2];
rz(2.5261836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0630546) q[1];
sx q[1];
rz(-1.6018493) q[1];
sx q[1];
rz(-1.6632776) q[1];
x q[2];
rz(-1.9014992) q[3];
sx q[3];
rz(-2.2281149) q[3];
sx q[3];
rz(0.98188321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25980514) q[2];
sx q[2];
rz(-1.2299579) q[2];
sx q[2];
rz(-1.7477431) q[2];
rz(-1.7143837) q[3];
sx q[3];
rz(-2.2673159) q[3];
sx q[3];
rz(-2.8425596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3706936) q[0];
sx q[0];
rz(-2.735266) q[0];
sx q[0];
rz(-2.4942177) q[0];
rz(2.8992843) q[1];
sx q[1];
rz(-1.4360042) q[1];
sx q[1];
rz(1.4586331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3577135) q[0];
sx q[0];
rz(-1.4365968) q[0];
sx q[0];
rz(2.8154897) q[0];
x q[1];
rz(2.3056612) q[2];
sx q[2];
rz(-2.3337337) q[2];
sx q[2];
rz(0.69617803) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9414166) q[1];
sx q[1];
rz(-1.5343018) q[1];
sx q[1];
rz(-1.0851651) q[1];
rz(-pi) q[2];
rz(0.69497739) q[3];
sx q[3];
rz(-0.91013346) q[3];
sx q[3];
rz(-1.7373178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4316537) q[2];
sx q[2];
rz(-0.84448758) q[2];
sx q[2];
rz(-1.9169774) q[2];
rz(0.25105181) q[3];
sx q[3];
rz(-0.71728388) q[3];
sx q[3];
rz(2.203598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5421903) q[0];
sx q[0];
rz(-2.0090065) q[0];
sx q[0];
rz(-0.19790025) q[0];
rz(-1.9056162) q[1];
sx q[1];
rz(-2.7695152) q[1];
sx q[1];
rz(3.1370251) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64070797) q[0];
sx q[0];
rz(-1.3237778) q[0];
sx q[0];
rz(2.5480854) q[0];
rz(-pi) q[1];
rz(1.0367583) q[2];
sx q[2];
rz(-1.6421841) q[2];
sx q[2];
rz(2.57881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98036042) q[1];
sx q[1];
rz(-1.2975818) q[1];
sx q[1];
rz(-1.4153764) q[1];
rz(-2.7907333) q[3];
sx q[3];
rz(-1.3727194) q[3];
sx q[3];
rz(2.4048865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68134754) q[2];
sx q[2];
rz(-0.494151) q[2];
sx q[2];
rz(2.8013308) q[2];
rz(-1.6081238) q[3];
sx q[3];
rz(-1.7483277) q[3];
sx q[3];
rz(1.5197915) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22572868) q[0];
sx q[0];
rz(-0.69750834) q[0];
sx q[0];
rz(-1.1874143) q[0];
rz(-1.7200708) q[1];
sx q[1];
rz(-2.7233796) q[1];
sx q[1];
rz(-3.0112867) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9604089) q[0];
sx q[0];
rz(-1.4474004) q[0];
sx q[0];
rz(0.058870319) q[0];
x q[1];
rz(-2.2557507) q[2];
sx q[2];
rz(-2.3315213) q[2];
sx q[2];
rz(2.5145234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7769717) q[1];
sx q[1];
rz(-1.4029619) q[1];
sx q[1];
rz(1.5701152) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6076419) q[3];
sx q[3];
rz(-0.26290694) q[3];
sx q[3];
rz(-2.2927566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80798951) q[2];
sx q[2];
rz(-1.5520059) q[2];
sx q[2];
rz(-1.0398593) q[2];
rz(-2.9676159) q[3];
sx q[3];
rz(-2.0128553) q[3];
sx q[3];
rz(-2.4719293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60798764) q[0];
sx q[0];
rz(-1.0796115) q[0];
sx q[0];
rz(0.84130353) q[0];
rz(-1.816642) q[1];
sx q[1];
rz(-2.1957896) q[1];
sx q[1];
rz(-1.673505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45638915) q[0];
sx q[0];
rz(-0.51420553) q[0];
sx q[0];
rz(-0.77001621) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1300342) q[2];
sx q[2];
rz(-1.0521787) q[2];
sx q[2];
rz(2.8457038) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74875301) q[1];
sx q[1];
rz(-2.3733099) q[1];
sx q[1];
rz(-0.70559754) q[1];
rz(-pi) q[2];
rz(1.2135336) q[3];
sx q[3];
rz(-1.8566086) q[3];
sx q[3];
rz(0.88102007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2738721) q[2];
sx q[2];
rz(-0.90110675) q[2];
sx q[2];
rz(-0.64485288) q[2];
rz(-1.0817179) q[3];
sx q[3];
rz(-2.7223301) q[3];
sx q[3];
rz(-0.7209512) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69407392) q[0];
sx q[0];
rz(-2.9677291) q[0];
sx q[0];
rz(-0.41282594) q[0];
rz(-1.5326356) q[1];
sx q[1];
rz(-0.91333476) q[1];
sx q[1];
rz(-2.434381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0577257) q[0];
sx q[0];
rz(-0.8208771) q[0];
sx q[0];
rz(0.42695257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4143432) q[2];
sx q[2];
rz(-2.3067637) q[2];
sx q[2];
rz(0.1231269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.096371264) q[1];
sx q[1];
rz(-2.1003083) q[1];
sx q[1];
rz(-1.96223) q[1];
rz(2.1096346) q[3];
sx q[3];
rz(-1.4488359) q[3];
sx q[3];
rz(1.7820047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1327208) q[2];
sx q[2];
rz(-2.8421695) q[2];
sx q[2];
rz(-0.51187619) q[2];
rz(-2.4608608) q[3];
sx q[3];
rz(-1.778089) q[3];
sx q[3];
rz(-0.16630047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7452288) q[0];
sx q[0];
rz(-1.5487211) q[0];
sx q[0];
rz(-2.6863099) q[0];
rz(1.0633172) q[1];
sx q[1];
rz(-0.89439193) q[1];
sx q[1];
rz(-3.0696226) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6980879) q[0];
sx q[0];
rz(-3.1060954) q[0];
sx q[0];
rz(0.3548236) q[0];
x q[1];
rz(-0.049111185) q[2];
sx q[2];
rz(-2.7002044) q[2];
sx q[2];
rz(0.15234824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1942206) q[1];
sx q[1];
rz(-0.76057077) q[1];
sx q[1];
rz(-1.0494227) q[1];
x q[2];
rz(2.2293363) q[3];
sx q[3];
rz(-2.7317762) q[3];
sx q[3];
rz(1.2959105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3528184) q[2];
sx q[2];
rz(-0.41663751) q[2];
sx q[2];
rz(-2.519506) q[2];
rz(-2.3173053) q[3];
sx q[3];
rz(-2.0485179) q[3];
sx q[3];
rz(-0.85339671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033009919) q[0];
sx q[0];
rz(-1.1065296) q[0];
sx q[0];
rz(-2.1906817) q[0];
rz(1.7225601) q[1];
sx q[1];
rz(-0.483069) q[1];
sx q[1];
rz(-0.61666617) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66835082) q[0];
sx q[0];
rz(-2.5010273) q[0];
sx q[0];
rz(0.040706445) q[0];
x q[1];
rz(1.1918814) q[2];
sx q[2];
rz(-0.99493631) q[2];
sx q[2];
rz(2.4873231) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2211799) q[1];
sx q[1];
rz(-1.3592459) q[1];
sx q[1];
rz(1.5736363) q[1];
x q[2];
rz(2.4343726) q[3];
sx q[3];
rz(-1.116718) q[3];
sx q[3];
rz(1.8994768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0509433) q[2];
sx q[2];
rz(-1.1530777) q[2];
sx q[2];
rz(-2.0743745) q[2];
rz(2.0459335) q[3];
sx q[3];
rz(-1.9039543) q[3];
sx q[3];
rz(0.9637951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.511956) q[0];
sx q[0];
rz(-2.1558599) q[0];
sx q[0];
rz(2.070367) q[0];
rz(1.4818954) q[1];
sx q[1];
rz(-1.0922468) q[1];
sx q[1];
rz(-2.560871) q[1];
rz(-0.26255519) q[2];
sx q[2];
rz(-1.5468183) q[2];
sx q[2];
rz(0.47894947) q[2];
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
