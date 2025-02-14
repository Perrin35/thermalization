OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(2.1535518) q[0];
sx q[0];
rz(9.0721985) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-3.92634) q[1];
sx q[1];
rz(9.3974455) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3685212) q[0];
sx q[0];
rz(-1.3734682) q[0];
sx q[0];
rz(-0.34223603) q[0];
rz(-1.4352082) q[2];
sx q[2];
rz(-1.7426881) q[2];
sx q[2];
rz(-2.7219049) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.084301) q[1];
sx q[1];
rz(-0.17070577) q[1];
sx q[1];
rz(-0.18526669) q[1];
rz(-pi) q[2];
rz(2.5952847) q[3];
sx q[3];
rz(-0.27537333) q[3];
sx q[3];
rz(-2.7978123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0693543) q[2];
sx q[2];
rz(-1.3212997) q[2];
sx q[2];
rz(1.5551152) q[2];
rz(0.72748264) q[3];
sx q[3];
rz(-0.95557094) q[3];
sx q[3];
rz(-0.47310841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6721866) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-2.8166215) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(2.6281338) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531813) q[0];
sx q[0];
rz(-1.9741749) q[0];
sx q[0];
rz(2.2405586) q[0];
x q[1];
rz(2.8186501) q[2];
sx q[2];
rz(-2.2720798) q[2];
sx q[2];
rz(0.27837929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70043515) q[1];
sx q[1];
rz(-2.0701906) q[1];
sx q[1];
rz(2.5446289) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88708892) q[3];
sx q[3];
rz(-2.4389431) q[3];
sx q[3];
rz(1.6446554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1309588) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(2.0531674) q[2];
rz(0.66793495) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(-1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208991) q[0];
sx q[0];
rz(-2.2405393) q[0];
sx q[0];
rz(1.2574842) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-3.0501922) q[1];
sx q[1];
rz(-0.65748293) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5705868) q[0];
sx q[0];
rz(-2.0181542) q[0];
sx q[0];
rz(-1.7156832) q[0];
x q[1];
rz(-2.3917603) q[2];
sx q[2];
rz(-0.74374357) q[2];
sx q[2];
rz(-2.2735571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0850683) q[1];
sx q[1];
rz(-2.6303362) q[1];
sx q[1];
rz(0.35447094) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90963962) q[3];
sx q[3];
rz(-1.3245021) q[3];
sx q[3];
rz(-1.1681739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(-1.8910889) q[2];
rz(1.6999543) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(2.3120248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.68469754) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(-2.5392505) q[0];
rz(0.96386987) q[1];
sx q[1];
rz(-2.3606221) q[1];
sx q[1];
rz(0.22183713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9070061) q[0];
sx q[0];
rz(-0.013628634) q[0];
sx q[0];
rz(2.5383294) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8922563) q[2];
sx q[2];
rz(-2.0866924) q[2];
sx q[2];
rz(0.14125401) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79961046) q[1];
sx q[1];
rz(-1.0162119) q[1];
sx q[1];
rz(-2.0444524) q[1];
rz(2.9939566) q[3];
sx q[3];
rz(-2.6872928) q[3];
sx q[3];
rz(2.4678178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43231371) q[2];
sx q[2];
rz(-1.1608492) q[2];
sx q[2];
rz(3.0529037) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(-3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13224193) q[0];
sx q[0];
rz(-0.049533822) q[0];
sx q[0];
rz(-0.90231878) q[0];
rz(-0.29014507) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(1.9754999) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571044) q[0];
sx q[0];
rz(-2.1323626) q[0];
sx q[0];
rz(1.7620371) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8140712) q[2];
sx q[2];
rz(-1.8406879) q[2];
sx q[2];
rz(-0.44651595) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8836959) q[1];
sx q[1];
rz(-1.9120815) q[1];
sx q[1];
rz(2.8410556) q[1];
x q[2];
rz(3.0405255) q[3];
sx q[3];
rz(-0.25288452) q[3];
sx q[3];
rz(-0.28836461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.12639283) q[2];
sx q[2];
rz(-0.79605278) q[2];
sx q[2];
rz(1.2882721) q[2];
rz(-2.1081693) q[3];
sx q[3];
rz(-2.0234334) q[3];
sx q[3];
rz(-1.6667295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7159336) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(-2.9820251) q[0];
rz(0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(-1.0608231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55218766) q[0];
sx q[0];
rz(-2.1656143) q[0];
sx q[0];
rz(-2.908559) q[0];
x q[1];
rz(2.873704) q[2];
sx q[2];
rz(-1.7888336) q[2];
sx q[2];
rz(1.6074558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6276803) q[1];
sx q[1];
rz(-0.6100756) q[1];
sx q[1];
rz(1.3327636) q[1];
rz(-1.4727598) q[3];
sx q[3];
rz(-2.1475002) q[3];
sx q[3];
rz(2.2876379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0173505) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(0.38786495) q[2];
rz(1.3346437) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(-2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(2.2482596) q[0];
rz(-0.54126254) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(1.062324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996595) q[0];
sx q[0];
rz(-2.3627653) q[0];
sx q[0];
rz(-2.7657484) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5028033) q[2];
sx q[2];
rz(-0.74882245) q[2];
sx q[2];
rz(-2.6523726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17328611) q[1];
sx q[1];
rz(-1.5943375) q[1];
sx q[1];
rz(0.55577718) q[1];
x q[2];
rz(-2.331953) q[3];
sx q[3];
rz(-1.7160549) q[3];
sx q[3];
rz(-0.64034772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58747753) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(-1.1678196) q[2];
rz(-2.3217412) q[3];
sx q[3];
rz(-0.76698118) q[3];
sx q[3];
rz(-0.64275536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.576936) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(-1.3078974) q[0];
rz(1.4564184) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(3.0110722) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042386656) q[0];
sx q[0];
rz(-1.4975647) q[0];
sx q[0];
rz(-1.8797148) q[0];
rz(-pi) q[1];
rz(-2.3442535) q[2];
sx q[2];
rz(-1.9249462) q[2];
sx q[2];
rz(2.7018099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3106892) q[1];
sx q[1];
rz(-0.92055087) q[1];
sx q[1];
rz(-1.6199153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7217312) q[3];
sx q[3];
rz(-0.85866195) q[3];
sx q[3];
rz(2.9809088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-0.83883494) q[2];
sx q[2];
rz(1.4107417) q[2];
rz(0.12061067) q[3];
sx q[3];
rz(-1.4863622) q[3];
sx q[3];
rz(1.6284778) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738344) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(3.0009785) q[0];
rz(2.4070542) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(-2.3908884) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0854617) q[0];
sx q[0];
rz(-1.3358677) q[0];
sx q[0];
rz(-0.30130193) q[0];
rz(-pi) q[1];
rz(-2.9042013) q[2];
sx q[2];
rz(-1.2331881) q[2];
sx q[2];
rz(0.92729502) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0486601) q[1];
sx q[1];
rz(-1.5092998) q[1];
sx q[1];
rz(1.7503121) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0680411) q[3];
sx q[3];
rz(-2.2760609) q[3];
sx q[3];
rz(-1.2770231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.004403) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(2.3889551) q[2];
rz(-0.33958069) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(-3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362157) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(-2.5035653) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(2.3400838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214695) q[0];
sx q[0];
rz(-0.82372181) q[0];
sx q[0];
rz(2.864564) q[0];
rz(-pi) q[1];
rz(1.837977) q[2];
sx q[2];
rz(-0.75088402) q[2];
sx q[2];
rz(0.24611404) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.90097016) q[1];
sx q[1];
rz(-2.3055446) q[1];
sx q[1];
rz(-0.052608629) q[1];
x q[2];
rz(-1.3574886) q[3];
sx q[3];
rz(-0.65591988) q[3];
sx q[3];
rz(1.3154958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8375497) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(-1.1207646) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-0.46375436) q[3];
sx q[3];
rz(-1.9161061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365702) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(0.90514056) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(-1.1726562) q[2];
sx q[2];
rz(-0.7672933) q[2];
sx q[2];
rz(2.2258022) q[2];
rz(1.0488679) q[3];
sx q[3];
rz(-1.3216139) q[3];
sx q[3];
rz(-0.25861964) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
