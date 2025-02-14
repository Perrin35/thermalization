OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0394734) q[0];
sx q[0];
rz(-1.4730299) q[0];
sx q[0];
rz(-3.0102475) q[0];
rz(0.036845358) q[1];
sx q[1];
rz(-0.52739066) q[1];
sx q[1];
rz(-1.3091458) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15802424) q[0];
sx q[0];
rz(-0.82497549) q[0];
sx q[0];
rz(1.0452861) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1151343) q[2];
sx q[2];
rz(-1.1725468) q[2];
sx q[2];
rz(-1.8894701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19840163) q[1];
sx q[1];
rz(-1.8089589) q[1];
sx q[1];
rz(-2.8203301) q[1];
x q[2];
rz(-2.3917624) q[3];
sx q[3];
rz(-0.24838003) q[3];
sx q[3];
rz(-1.6102143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5347791) q[2];
sx q[2];
rz(-2.7853192) q[2];
sx q[2];
rz(2.479539) q[2];
rz(0.74696294) q[3];
sx q[3];
rz(-2.2253939) q[3];
sx q[3];
rz(1.0158739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5544283) q[0];
sx q[0];
rz(-2.9968408) q[0];
sx q[0];
rz(-1.9594877) q[0];
rz(-2.3960522) q[1];
sx q[1];
rz(-2.5423971) q[1];
sx q[1];
rz(-3.0381957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1348159) q[0];
sx q[0];
rz(-2.6511483) q[0];
sx q[0];
rz(-2.6314526) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1680122) q[2];
sx q[2];
rz(-1.2478618) q[2];
sx q[2];
rz(-2.3514049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34161257) q[1];
sx q[1];
rz(-0.80474058) q[1];
sx q[1];
rz(-1.1031723) q[1];
rz(-pi) q[2];
rz(1.7595791) q[3];
sx q[3];
rz(-0.75136614) q[3];
sx q[3];
rz(2.7392859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5891002) q[2];
sx q[2];
rz(-1.8742259) q[2];
sx q[2];
rz(2.6972771) q[2];
rz(1.0537423) q[3];
sx q[3];
rz(-0.070662347) q[3];
sx q[3];
rz(1.1963371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0788197) q[0];
sx q[0];
rz(-1.2508996) q[0];
sx q[0];
rz(0.66251063) q[0];
rz(-1.6569116) q[1];
sx q[1];
rz(-2.1187481) q[1];
sx q[1];
rz(2.9248765) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92531918) q[0];
sx q[0];
rz(-2.0331622) q[0];
sx q[0];
rz(-2.3564649) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6609801) q[2];
sx q[2];
rz(-0.79250249) q[2];
sx q[2];
rz(-2.8414244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3308476) q[1];
sx q[1];
rz(-2.243089) q[1];
sx q[1];
rz(-3.0903893) q[1];
x q[2];
rz(-1.7314579) q[3];
sx q[3];
rz(-1.0815797) q[3];
sx q[3];
rz(0.9582203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66389877) q[2];
sx q[2];
rz(-2.615216) q[2];
sx q[2];
rz(2.9065175) q[2];
rz(-2.7663686) q[3];
sx q[3];
rz(-0.90685654) q[3];
sx q[3];
rz(-0.71845976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5821871) q[0];
sx q[0];
rz(-0.087787293) q[0];
sx q[0];
rz(-1.0002332) q[0];
rz(-1.2031215) q[1];
sx q[1];
rz(-1.5088046) q[1];
sx q[1];
rz(2.5968754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398374) q[0];
sx q[0];
rz(-2.4271963) q[0];
sx q[0];
rz(-0.06846662) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3339231) q[2];
sx q[2];
rz(-0.90014202) q[2];
sx q[2];
rz(-2.9602221) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96918584) q[1];
sx q[1];
rz(-1.6891857) q[1];
sx q[1];
rz(-1.7612996) q[1];
rz(-0.95787134) q[3];
sx q[3];
rz(-1.8184693) q[3];
sx q[3];
rz(-1.6771468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41877052) q[2];
sx q[2];
rz(-0.79221574) q[2];
sx q[2];
rz(0.79696068) q[2];
rz(-0.289251) q[3];
sx q[3];
rz(-2.1713493) q[3];
sx q[3];
rz(-2.7204035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5274984) q[0];
sx q[0];
rz(-3.0526563) q[0];
sx q[0];
rz(1.2689137) q[0];
rz(0.96619636) q[1];
sx q[1];
rz(-1.7485917) q[1];
sx q[1];
rz(-0.46636811) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0206576) q[0];
sx q[0];
rz(-1.0758721) q[0];
sx q[0];
rz(-1.0878956) q[0];
rz(-pi) q[1];
rz(0.46196533) q[2];
sx q[2];
rz(-1.8868539) q[2];
sx q[2];
rz(2.8082928) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5261698) q[1];
sx q[1];
rz(-0.96994441) q[1];
sx q[1];
rz(2.7533965) q[1];
x q[2];
rz(-0.39042274) q[3];
sx q[3];
rz(-2.1964034) q[3];
sx q[3];
rz(2.1805891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40631488) q[2];
sx q[2];
rz(-2.3390528) q[2];
sx q[2];
rz(2.3373513) q[2];
rz(-2.6546226) q[3];
sx q[3];
rz(-2.099497) q[3];
sx q[3];
rz(-3.13412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93953472) q[0];
sx q[0];
rz(-0.29629961) q[0];
sx q[0];
rz(-0.95788389) q[0];
rz(1.6288039) q[1];
sx q[1];
rz(-2.084338) q[1];
sx q[1];
rz(-1.588795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4544425) q[0];
sx q[0];
rz(-1.5936562) q[0];
sx q[0];
rz(-1.6506565) q[0];
x q[1];
rz(-0.79926305) q[2];
sx q[2];
rz(-2.7529319) q[2];
sx q[2];
rz(-0.51539076) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54974906) q[1];
sx q[1];
rz(-1.3938245) q[1];
sx q[1];
rz(-0.66284499) q[1];
x q[2];
rz(-0.62726043) q[3];
sx q[3];
rz(-1.2377487) q[3];
sx q[3];
rz(1.8740538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2019041) q[2];
sx q[2];
rz(-1.0558015) q[2];
sx q[2];
rz(-0.73516694) q[2];
rz(-0.9489263) q[3];
sx q[3];
rz(-0.85845033) q[3];
sx q[3];
rz(0.86400509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76148024) q[0];
sx q[0];
rz(-1.004847) q[0];
sx q[0];
rz(2.5349706) q[0];
rz(2.3573549) q[1];
sx q[1];
rz(-1.3332858) q[1];
sx q[1];
rz(-1.3444791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30596581) q[0];
sx q[0];
rz(-2.2971662) q[0];
sx q[0];
rz(-2.4854922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6450301) q[2];
sx q[2];
rz(-1.3482058) q[2];
sx q[2];
rz(-0.45997657) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.231338) q[1];
sx q[1];
rz(-1.7926072) q[1];
sx q[1];
rz(0.40738687) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.090769692) q[3];
sx q[3];
rz(-1.1072973) q[3];
sx q[3];
rz(0.32583729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.532423) q[2];
sx q[2];
rz(-0.50749856) q[2];
sx q[2];
rz(-1.3896821) q[2];
rz(-0.17942795) q[3];
sx q[3];
rz(-2.1556985) q[3];
sx q[3];
rz(-0.40670407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0591902) q[0];
sx q[0];
rz(-0.32408369) q[0];
sx q[0];
rz(2.3865336) q[0];
rz(-0.45626196) q[1];
sx q[1];
rz(-1.4249964) q[1];
sx q[1];
rz(2.0645352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2307325) q[0];
sx q[0];
rz(-2.5344779) q[0];
sx q[0];
rz(-0.40672983) q[0];
rz(-pi) q[1];
rz(-1.6918534) q[2];
sx q[2];
rz(-0.79059764) q[2];
sx q[2];
rz(-1.3644753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7491319) q[1];
sx q[1];
rz(-1.2866255) q[1];
sx q[1];
rz(-2.6662213) q[1];
rz(-2.9829426) q[3];
sx q[3];
rz(-2.1057743) q[3];
sx q[3];
rz(-1.3033997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72851744) q[2];
sx q[2];
rz(-2.04144) q[2];
sx q[2];
rz(2.4033578) q[2];
rz(-1.5089367) q[3];
sx q[3];
rz(-1.3239219) q[3];
sx q[3];
rz(-1.1608605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64410925) q[0];
sx q[0];
rz(-1.4583541) q[0];
sx q[0];
rz(2.4926376) q[0];
rz(0.68967462) q[1];
sx q[1];
rz(-0.70976218) q[1];
sx q[1];
rz(-2.7241657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6605292) q[0];
sx q[0];
rz(-1.3205577) q[0];
sx q[0];
rz(2.5323243) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8418188) q[2];
sx q[2];
rz(-1.4978786) q[2];
sx q[2];
rz(0.48540965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46681225) q[1];
sx q[1];
rz(-1.1622475) q[1];
sx q[1];
rz(-1.684434) q[1];
rz(0.77885742) q[3];
sx q[3];
rz(-0.53611272) q[3];
sx q[3];
rz(-0.4975618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6680341) q[2];
sx q[2];
rz(-1.4213976) q[2];
sx q[2];
rz(-3.0958946) q[2];
rz(2.8081196) q[3];
sx q[3];
rz(-2.5662751) q[3];
sx q[3];
rz(0.12588178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4561653) q[0];
sx q[0];
rz(-0.5394772) q[0];
sx q[0];
rz(1.217655) q[0];
rz(-0.24670163) q[1];
sx q[1];
rz(-1.291899) q[1];
sx q[1];
rz(0.47952476) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8586774) q[0];
sx q[0];
rz(-2.0020131) q[0];
sx q[0];
rz(-2.2515576) q[0];
rz(-pi) q[1];
rz(0.99202435) q[2];
sx q[2];
rz(-2.357748) q[2];
sx q[2];
rz(-1.0933361) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.785276) q[1];
sx q[1];
rz(-1.061376) q[1];
sx q[1];
rz(1.9059577) q[1];
rz(-pi) q[2];
rz(0.44879313) q[3];
sx q[3];
rz(-2.6330607) q[3];
sx q[3];
rz(-0.75565805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55629998) q[2];
sx q[2];
rz(-1.2703398) q[2];
sx q[2];
rz(-1.7362107) q[2];
rz(2.4893238) q[3];
sx q[3];
rz(-2.8758719) q[3];
sx q[3];
rz(-0.71776596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2300867) q[0];
sx q[0];
rz(-1.7510887) q[0];
sx q[0];
rz(-0.44816309) q[0];
rz(0.74408342) q[1];
sx q[1];
rz(-2.4241445) q[1];
sx q[1];
rz(-1.4345899) q[1];
rz(-2.9058331) q[2];
sx q[2];
rz(-2.5452062) q[2];
sx q[2];
rz(1.7231981) q[2];
rz(-2.1461829) q[3];
sx q[3];
rz(-1.2708455) q[3];
sx q[3];
rz(1.3686913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
