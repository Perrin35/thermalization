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
rz(0.67233664) q[0];
sx q[0];
rz(-1.2824751) q[0];
sx q[0];
rz(1.6096492) q[0];
rz(1.2338282) q[1];
sx q[1];
rz(-1.9246074) q[1];
sx q[1];
rz(-1.4438862) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1506526) q[0];
sx q[0];
rz(-0.67545891) q[0];
sx q[0];
rz(1.4399685) q[0];
rz(0.89272372) q[2];
sx q[2];
rz(-2.6916382) q[2];
sx q[2];
rz(-1.8261248) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.68221625) q[1];
sx q[1];
rz(-2.3704048) q[1];
sx q[1];
rz(-2.5384063) q[1];
rz(-1.6495632) q[3];
sx q[3];
rz(-0.34011671) q[3];
sx q[3];
rz(-0.3398557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3621138) q[2];
sx q[2];
rz(-1.9356091) q[2];
sx q[2];
rz(-0.77902478) q[2];
rz(1.834076) q[3];
sx q[3];
rz(-0.3236168) q[3];
sx q[3];
rz(-1.754508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5883412) q[0];
sx q[0];
rz(-1.7760176) q[0];
sx q[0];
rz(-0.27728444) q[0];
rz(-0.39018997) q[1];
sx q[1];
rz(-2.0398102) q[1];
sx q[1];
rz(2.0139093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4235315) q[0];
sx q[0];
rz(-0.047022659) q[0];
sx q[0];
rz(0.440098) q[0];
x q[1];
rz(0.28553005) q[2];
sx q[2];
rz(-1.6678311) q[2];
sx q[2];
rz(-0.75733987) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8133328) q[1];
sx q[1];
rz(-1.5082803) q[1];
sx q[1];
rz(-1.2986138) q[1];
x q[2];
rz(0.66172285) q[3];
sx q[3];
rz(-1.265101) q[3];
sx q[3];
rz(0.88246417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7057425) q[2];
sx q[2];
rz(-2.5472992) q[2];
sx q[2];
rz(-1.1951813) q[2];
rz(-0.65520203) q[3];
sx q[3];
rz(-1.4846669) q[3];
sx q[3];
rz(-0.5932194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6708267) q[0];
sx q[0];
rz(-2.4556181) q[0];
sx q[0];
rz(0.92054787) q[0];
rz(1.8596733) q[1];
sx q[1];
rz(-2.0125407) q[1];
sx q[1];
rz(-1.4528073) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1143415) q[0];
sx q[0];
rz(-0.14273164) q[0];
sx q[0];
rz(0.65988512) q[0];
rz(-0.31450504) q[2];
sx q[2];
rz(-2.1172519) q[2];
sx q[2];
rz(0.98775452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5389899) q[1];
sx q[1];
rz(-0.33947152) q[1];
sx q[1];
rz(2.9871639) q[1];
rz(-0.7887855) q[3];
sx q[3];
rz(-0.093063912) q[3];
sx q[3];
rz(-2.399202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.526223) q[2];
sx q[2];
rz(-1.0542032) q[2];
sx q[2];
rz(2.7426381) q[2];
rz(-2.6279972) q[3];
sx q[3];
rz(-0.31702888) q[3];
sx q[3];
rz(1.1494466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32617635) q[0];
sx q[0];
rz(-0.57891095) q[0];
sx q[0];
rz(-1.2166566) q[0];
rz(-0.51219621) q[1];
sx q[1];
rz(-1.1108111) q[1];
sx q[1];
rz(-1.1042575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35262689) q[0];
sx q[0];
rz(-1.2944844) q[0];
sx q[0];
rz(0.31310149) q[0];
rz(-pi) q[1];
rz(-2.9680164) q[2];
sx q[2];
rz(-1.3921228) q[2];
sx q[2];
rz(1.0965958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.11286375) q[1];
sx q[1];
rz(-0.20564889) q[1];
sx q[1];
rz(1.8702694) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8301303) q[3];
sx q[3];
rz(-1.7899198) q[3];
sx q[3];
rz(0.27838368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5451374) q[2];
sx q[2];
rz(-1.1918273) q[2];
sx q[2];
rz(2.5917501) q[2];
rz(2.1846629) q[3];
sx q[3];
rz(-0.79920971) q[3];
sx q[3];
rz(1.6526615) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4855708) q[0];
sx q[0];
rz(-2.2915884) q[0];
sx q[0];
rz(2.3824084) q[0];
rz(0.94866577) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(-2.8099828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054197) q[0];
sx q[0];
rz(-2.3080462) q[0];
sx q[0];
rz(1.1655318) q[0];
x q[1];
rz(-2.6668915) q[2];
sx q[2];
rz(-1.6031638) q[2];
sx q[2];
rz(-0.84454483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75237235) q[1];
sx q[1];
rz(-1.3942317) q[1];
sx q[1];
rz(0.15644381) q[1];
rz(-pi) q[2];
rz(1.4344352) q[3];
sx q[3];
rz(-1.162774) q[3];
sx q[3];
rz(-1.3566164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9552976) q[2];
sx q[2];
rz(-0.95185995) q[2];
sx q[2];
rz(-0.9430421) q[2];
rz(-2.9694563) q[3];
sx q[3];
rz(-2.19682) q[3];
sx q[3];
rz(0.22323639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.405769) q[0];
sx q[0];
rz(-0.74546927) q[0];
sx q[0];
rz(1.5752342) q[0];
rz(0.15815059) q[1];
sx q[1];
rz(-2.3731396) q[1];
sx q[1];
rz(-2.2134773) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.428663) q[0];
sx q[0];
rz(-0.33131772) q[0];
sx q[0];
rz(0.40508799) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42144962) q[2];
sx q[2];
rz(-1.0578276) q[2];
sx q[2];
rz(-2.9861272) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0911407) q[1];
sx q[1];
rz(-2.4311922) q[1];
sx q[1];
rz(-1.1732167) q[1];
rz(-1.6949878) q[3];
sx q[3];
rz(-1.2203794) q[3];
sx q[3];
rz(-2.7563376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1930262) q[2];
sx q[2];
rz(-2.8552449) q[2];
sx q[2];
rz(-2.7890653) q[2];
rz(-2.7726717) q[3];
sx q[3];
rz(-0.58266321) q[3];
sx q[3];
rz(-0.85744706) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7848659) q[0];
sx q[0];
rz(-3.1251188) q[0];
sx q[0];
rz(-2.4705868) q[0];
rz(-0.10637936) q[1];
sx q[1];
rz(-2.2439067) q[1];
sx q[1];
rz(2.5118714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2248005) q[0];
sx q[0];
rz(-2.9491701) q[0];
sx q[0];
rz(0.097956603) q[0];
x q[1];
rz(2.7313298) q[2];
sx q[2];
rz(-1.9681491) q[2];
sx q[2];
rz(0.58471459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56901819) q[1];
sx q[1];
rz(-0.33199939) q[1];
sx q[1];
rz(-0.23004679) q[1];
x q[2];
rz(-1.2604146) q[3];
sx q[3];
rz(-2.2058626) q[3];
sx q[3];
rz(0.47582754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.070039198) q[2];
sx q[2];
rz(-2.2420501) q[2];
sx q[2];
rz(-2.1001935) q[2];
rz(1.6085767) q[3];
sx q[3];
rz(-2.2074015) q[3];
sx q[3];
rz(1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94569412) q[0];
sx q[0];
rz(-0.65106374) q[0];
sx q[0];
rz(-2.661929) q[0];
rz(3.116963) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(0.13951313) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280645) q[0];
sx q[0];
rz(-1.6798927) q[0];
sx q[0];
rz(-0.086221545) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.191338) q[2];
sx q[2];
rz(-0.59796158) q[2];
sx q[2];
rz(0.68956748) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0932659) q[1];
sx q[1];
rz(-0.67950002) q[1];
sx q[1];
rz(-0.26643403) q[1];
rz(-3.0728389) q[3];
sx q[3];
rz(-2.4700751) q[3];
sx q[3];
rz(2.9686787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8861683) q[2];
sx q[2];
rz(-1.5443799) q[2];
sx q[2];
rz(-1.1523979) q[2];
rz(-1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(1.8324119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.060147978) q[0];
sx q[0];
rz(-1.6763433) q[0];
sx q[0];
rz(-1.4981221) q[0];
rz(-0.76088798) q[1];
sx q[1];
rz(-2.1983169) q[1];
sx q[1];
rz(-1.6852185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8955904) q[0];
sx q[0];
rz(-2.1293679) q[0];
sx q[0];
rz(-1.7825401) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4001021) q[2];
sx q[2];
rz(-1.7662373) q[2];
sx q[2];
rz(-0.2646499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0800164) q[1];
sx q[1];
rz(-1.9655088) q[1];
sx q[1];
rz(1.030974) q[1];
x q[2];
rz(2.7237503) q[3];
sx q[3];
rz(-2.0977694) q[3];
sx q[3];
rz(-2.9067519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3219519) q[2];
sx q[2];
rz(-2.7072622) q[2];
sx q[2];
rz(-1.0158018) q[2];
rz(-0.89456931) q[3];
sx q[3];
rz(-1.2457448) q[3];
sx q[3];
rz(-2.3999124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9858518) q[0];
sx q[0];
rz(-2.4116801) q[0];
sx q[0];
rz(2.1515382) q[0];
rz(0.92297018) q[1];
sx q[1];
rz(-1.3786517) q[1];
sx q[1];
rz(-0.96614456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21912757) q[0];
sx q[0];
rz(-0.12258633) q[0];
sx q[0];
rz(-2.1845093) q[0];
x q[1];
rz(-0.43010148) q[2];
sx q[2];
rz(-1.9942339) q[2];
sx q[2];
rz(0.11030876) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50033874) q[1];
sx q[1];
rz(-0.49420824) q[1];
sx q[1];
rz(-0.92730913) q[1];
x q[2];
rz(2.0851785) q[3];
sx q[3];
rz(-0.22795151) q[3];
sx q[3];
rz(-2.0755656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82603377) q[2];
sx q[2];
rz(-2.3580599) q[2];
sx q[2];
rz(0.75251904) q[2];
rz(0.13713947) q[3];
sx q[3];
rz(-0.64118853) q[3];
sx q[3];
rz(-2.7916059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8990477) q[0];
sx q[0];
rz(-0.62472961) q[0];
sx q[0];
rz(-0.20150264) q[0];
rz(1.7789727) q[1];
sx q[1];
rz(-0.87281223) q[1];
sx q[1];
rz(2.9760117) q[1];
rz(-2.0803484) q[2];
sx q[2];
rz(-3.0399418) q[2];
sx q[2];
rz(-2.1352482) q[2];
rz(-1.9769382) q[3];
sx q[3];
rz(-1.9057169) q[3];
sx q[3];
rz(2.9939772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
