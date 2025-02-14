OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4234023) q[0];
sx q[0];
rz(-0.0027522491) q[0];
sx q[0];
rz(2.1299074) q[0];
rz(0.45971316) q[1];
sx q[1];
rz(4.4432321) q[1];
sx q[1];
rz(9.2530773) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3637488) q[0];
sx q[0];
rz(-2.1743589) q[0];
sx q[0];
rz(1.4715172) q[0];
rz(0.5131716) q[2];
sx q[2];
rz(-2.9193239) q[2];
sx q[2];
rz(-1.5602416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5682809) q[1];
sx q[1];
rz(-2.3900095) q[1];
sx q[1];
rz(0.64935301) q[1];
rz(-1.6640673) q[3];
sx q[3];
rz(-1.516946) q[3];
sx q[3];
rz(-0.71291956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1823938) q[2];
sx q[2];
rz(-0.48754075) q[2];
sx q[2];
rz(1.4200776) q[2];
rz(-1.9567664) q[3];
sx q[3];
rz(-1.5337557) q[3];
sx q[3];
rz(-2.0737295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0576393) q[0];
sx q[0];
rz(-1.6211442) q[0];
sx q[0];
rz(-2.6481096) q[0];
rz(-2.3656942) q[1];
sx q[1];
rz(-0.50965613) q[1];
sx q[1];
rz(0.50481558) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1238839) q[0];
sx q[0];
rz(-1.0320745) q[0];
sx q[0];
rz(-2.3305642) q[0];
rz(-pi) q[1];
rz(1.8680598) q[2];
sx q[2];
rz(-1.0631764) q[2];
sx q[2];
rz(2.076869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3146554) q[1];
sx q[1];
rz(-2.3874904) q[1];
sx q[1];
rz(0.60390632) q[1];
rz(0.26594578) q[3];
sx q[3];
rz(-1.5107181) q[3];
sx q[3];
rz(1.657965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33727553) q[2];
sx q[2];
rz(-1.8307468) q[2];
sx q[2];
rz(2.6709225) q[2];
rz(-2.5110631) q[3];
sx q[3];
rz(-1.0258521) q[3];
sx q[3];
rz(-1.4001018) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640963) q[0];
sx q[0];
rz(-0.12501669) q[0];
sx q[0];
rz(2.4858544) q[0];
rz(0.31133044) q[1];
sx q[1];
rz(-2.2475524) q[1];
sx q[1];
rz(-0.071455926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8250803) q[0];
sx q[0];
rz(-1.5208172) q[0];
sx q[0];
rz(-1.9216538) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3089789) q[2];
sx q[2];
rz(-0.48093167) q[2];
sx q[2];
rz(2.0455751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1972547) q[1];
sx q[1];
rz(-1.6363012) q[1];
sx q[1];
rz(1.08169) q[1];
rz(-pi) q[2];
rz(2.7903665) q[3];
sx q[3];
rz(-0.5595419) q[3];
sx q[3];
rz(2.9282687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91325703) q[2];
sx q[2];
rz(-1.5847289) q[2];
sx q[2];
rz(3.081591) q[2];
rz(-1.2126728) q[3];
sx q[3];
rz(-2.4433177) q[3];
sx q[3];
rz(2.3771299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54670984) q[0];
sx q[0];
rz(-1.8130274) q[0];
sx q[0];
rz(-2.5643964) q[0];
rz(-2.3120841) q[1];
sx q[1];
rz(-1.0521768) q[1];
sx q[1];
rz(-0.83782354) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5853607) q[0];
sx q[0];
rz(-0.62667003) q[0];
sx q[0];
rz(3.0100432) q[0];
rz(-0.060537593) q[2];
sx q[2];
rz(-2.1601387) q[2];
sx q[2];
rz(2.4174487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7260704) q[1];
sx q[1];
rz(-0.94968098) q[1];
sx q[1];
rz(0.46393053) q[1];
rz(2.1076848) q[3];
sx q[3];
rz(-2.3903525) q[3];
sx q[3];
rz(-0.43207016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.025658) q[2];
sx q[2];
rz(-1.5237153) q[2];
sx q[2];
rz(1.2238097) q[2];
rz(-0.4661679) q[3];
sx q[3];
rz(-2.7397082) q[3];
sx q[3];
rz(-2.536072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25662988) q[0];
sx q[0];
rz(-1.7054568) q[0];
sx q[0];
rz(-0.10109854) q[0];
rz(-2.7283607) q[1];
sx q[1];
rz(-2.7006472) q[1];
sx q[1];
rz(-2.0535927) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71095419) q[0];
sx q[0];
rz(-2.0530409) q[0];
sx q[0];
rz(2.7706465) q[0];
x q[1];
rz(1.1674983) q[2];
sx q[2];
rz(-2.3638862) q[2];
sx q[2];
rz(2.4571927) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6551825) q[1];
sx q[1];
rz(-1.2401199) q[1];
sx q[1];
rz(-2.8831456) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66941485) q[3];
sx q[3];
rz(-1.515404) q[3];
sx q[3];
rz(-2.9852569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92855144) q[2];
sx q[2];
rz(-2.2553359) q[2];
sx q[2];
rz(-1.586033) q[2];
rz(-1.0726311) q[3];
sx q[3];
rz(-1.6173247) q[3];
sx q[3];
rz(-0.13256375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9685386) q[0];
sx q[0];
rz(-0.88468495) q[0];
sx q[0];
rz(1.435085) q[0];
rz(0.75812078) q[1];
sx q[1];
rz(-0.70223141) q[1];
sx q[1];
rz(-3.039956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1124737) q[0];
sx q[0];
rz(-2.2879507) q[0];
sx q[0];
rz(2.4957335) q[0];
rz(-pi) q[1];
rz(0.91288699) q[2];
sx q[2];
rz(-0.55333558) q[2];
sx q[2];
rz(2.0408415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1730301) q[1];
sx q[1];
rz(-2.8391339) q[1];
sx q[1];
rz(-2.7773501) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6061344) q[3];
sx q[3];
rz(-1.0046008) q[3];
sx q[3];
rz(-2.9009852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8018084) q[2];
sx q[2];
rz(-2.0399751) q[2];
sx q[2];
rz(-1.3327117) q[2];
rz(1.0821651) q[3];
sx q[3];
rz(-2.3887631) q[3];
sx q[3];
rz(1.7180721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69671714) q[0];
sx q[0];
rz(-1.2514021) q[0];
sx q[0];
rz(-0.74309293) q[0];
rz(-2.3785059) q[1];
sx q[1];
rz(-2.5021195) q[1];
sx q[1];
rz(2.2230164) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7189499) q[0];
sx q[0];
rz(-0.026233679) q[0];
sx q[0];
rz(-1.2573186) q[0];
rz(-pi) q[1];
rz(-2.1821676) q[2];
sx q[2];
rz(-0.55771962) q[2];
sx q[2];
rz(0.096867933) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9715226) q[1];
sx q[1];
rz(-2.5405209) q[1];
sx q[1];
rz(-0.55723377) q[1];
x q[2];
rz(0.72414805) q[3];
sx q[3];
rz(-1.482943) q[3];
sx q[3];
rz(-1.4081362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4172198) q[2];
sx q[2];
rz(-1.1823187) q[2];
sx q[2];
rz(2.7057538) q[2];
rz(3.0271652) q[3];
sx q[3];
rz(-2.8947713) q[3];
sx q[3];
rz(0.59823263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.33335394) q[0];
sx q[0];
rz(-2.5840608) q[0];
sx q[0];
rz(0.58407855) q[0];
rz(-2.7059879) q[1];
sx q[1];
rz(-1.4608811) q[1];
sx q[1];
rz(-1.5199419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95660644) q[0];
sx q[0];
rz(-1.865956) q[0];
sx q[0];
rz(-1.0888238) q[0];
x q[1];
rz(0.21824117) q[2];
sx q[2];
rz(-1.1704966) q[2];
sx q[2];
rz(-1.9298273) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78970026) q[1];
sx q[1];
rz(-1.5613817) q[1];
sx q[1];
rz(-0.20603754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.734821) q[3];
sx q[3];
rz(-1.3151405) q[3];
sx q[3];
rz(-2.260316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0674151) q[2];
sx q[2];
rz(-1.7895074) q[2];
sx q[2];
rz(2.0254859) q[2];
rz(-1.3793147) q[3];
sx q[3];
rz(-2.0383056) q[3];
sx q[3];
rz(3.0232159) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922358) q[0];
sx q[0];
rz(-3.0638969) q[0];
sx q[0];
rz(2.351601) q[0];
rz(-3.0780011) q[1];
sx q[1];
rz(-0.60706943) q[1];
sx q[1];
rz(-2.7008609) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9544308) q[0];
sx q[0];
rz(-1.8858028) q[0];
sx q[0];
rz(-0.31401547) q[0];
rz(0.67040261) q[2];
sx q[2];
rz(-2.5385661) q[2];
sx q[2];
rz(-2.5489901) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1524618) q[1];
sx q[1];
rz(-0.50639443) q[1];
sx q[1];
rz(1.9146862) q[1];
rz(-pi) q[2];
rz(-1.7601898) q[3];
sx q[3];
rz(-1.6114858) q[3];
sx q[3];
rz(0.19408801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8792087) q[2];
sx q[2];
rz(-2.3449506) q[2];
sx q[2];
rz(0.66514307) q[2];
rz(-0.77196676) q[3];
sx q[3];
rz(-2.3699581) q[3];
sx q[3];
rz(1.5099248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7803698) q[0];
sx q[0];
rz(-2.4936115) q[0];
sx q[0];
rz(-1.004647) q[0];
rz(1.4077582) q[1];
sx q[1];
rz(-1.6561457) q[1];
sx q[1];
rz(1.8515324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697851) q[0];
sx q[0];
rz(-1.3243073) q[0];
sx q[0];
rz(2.3355961) q[0];
rz(-pi) q[1];
rz(0.97083135) q[2];
sx q[2];
rz(-1.7391676) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.066034868) q[1];
sx q[1];
rz(-1.3595288) q[1];
sx q[1];
rz(2.8610364) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96995207) q[3];
sx q[3];
rz(-1.4267216) q[3];
sx q[3];
rz(-0.010330095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1349692) q[2];
sx q[2];
rz(-1.9776191) q[2];
sx q[2];
rz(-0.25035614) q[2];
rz(-2.1641459) q[3];
sx q[3];
rz(-1.2075295) q[3];
sx q[3];
rz(-1.4704963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95494315) q[0];
sx q[0];
rz(-1.432812) q[0];
sx q[0];
rz(-1.1695255) q[0];
rz(-0.74465887) q[1];
sx q[1];
rz(-2.3458993) q[1];
sx q[1];
rz(0.69534272) q[1];
rz(1.5834687) q[2];
sx q[2];
rz(-1.7141533) q[2];
sx q[2];
rz(-1.2319596) q[2];
rz(2.70784) q[3];
sx q[3];
rz(-2.0278145) q[3];
sx q[3];
rz(-0.47587004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
