OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33485246) q[2];
sx q[2];
rz(-1.9960253) q[2];
sx q[2];
rz(1.877117) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4760121) q[1];
sx q[1];
rz(-1.952938) q[1];
sx q[1];
rz(1.0784472) q[1];
x q[2];
rz(0.11043926) q[3];
sx q[3];
rz(-0.38023284) q[3];
sx q[3];
rz(1.5649753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-0.28960323) q[2];
rz(2.2662207) q[3];
sx q[3];
rz(-1.0062199) q[3];
sx q[3];
rz(3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61525476) q[0];
sx q[0];
rz(-2.2741788) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.7864236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25027572) q[0];
sx q[0];
rz(-2.180897) q[0];
sx q[0];
rz(-0.079205714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6927035) q[2];
sx q[2];
rz(-1.7322455) q[2];
sx q[2];
rz(-0.09300692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.558074) q[1];
sx q[1];
rz(-2.0194224) q[1];
sx q[1];
rz(0.070986991) q[1];
rz(-pi) q[2];
rz(1.9676898) q[3];
sx q[3];
rz(-2.9429884) q[3];
sx q[3];
rz(1.393115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(0.53768349) q[2];
rz(-2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(2.1311549) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66353345) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-2.5007201) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(1.5887567) q[0];
x q[1];
rz(-2.2674019) q[2];
sx q[2];
rz(-2.0421931) q[2];
sx q[2];
rz(0.092560571) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6667337) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(-0.95400793) q[1];
rz(-2.8912656) q[3];
sx q[3];
rz(-2.2776789) q[3];
sx q[3];
rz(-1.7771429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(0.71737814) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(3.1304741) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3664368) q[0];
sx q[0];
rz(-1.1661134) q[0];
sx q[0];
rz(-1.0767656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83799329) q[2];
sx q[2];
rz(-2.0760771) q[2];
sx q[2];
rz(2.879564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30918446) q[1];
sx q[1];
rz(-2.5384181) q[1];
sx q[1];
rz(-2.1716154) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6075881) q[3];
sx q[3];
rz(-0.36557331) q[3];
sx q[3];
rz(0.064985736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.49542385) q[2];
sx q[2];
rz(-1.8858706) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-0.59316558) q[3];
sx q[3];
rz(-0.88808131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0304612) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(0.49452531) q[1];
sx q[1];
rz(-1.2978413) q[1];
sx q[1];
rz(-1.8146851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4715695) q[0];
sx q[0];
rz(-1.7243392) q[0];
sx q[0];
rz(-2.3877386) q[0];
rz(-pi) q[1];
rz(-0.58381501) q[2];
sx q[2];
rz(-0.8875672) q[2];
sx q[2];
rz(-0.77606397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(-2.521442) q[1];
x q[2];
rz(0.29019659) q[3];
sx q[3];
rz(-1.1786596) q[3];
sx q[3];
rz(-0.41662595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.999324) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.8959321) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-3.1384387) q[0];
sx q[0];
rz(-2.4601049) q[0];
rz(0.20755126) q[1];
sx q[1];
rz(-2.6782268) q[1];
sx q[1];
rz(-1.1157657) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(2.1372165) q[0];
rz(-pi) q[1];
rz(2.7115466) q[2];
sx q[2];
rz(-2.1950245) q[2];
sx q[2];
rz(2.4732694) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4286194) q[1];
sx q[1];
rz(-1.5884807) q[1];
sx q[1];
rz(-0.58847217) q[1];
rz(-1.7498155) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(-0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21268022) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(-0.35432717) q[2];
rz(0.3195233) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0091244) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-0.66450417) q[1];
sx q[1];
rz(-2.5792714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494941) q[0];
sx q[0];
rz(-1.5856165) q[0];
sx q[0];
rz(-0.90256079) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.855905) q[2];
sx q[2];
rz(-0.34491587) q[2];
sx q[2];
rz(0.069318511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6866236) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.7186233) q[1];
rz(-pi) q[2];
rz(1.4427156) q[3];
sx q[3];
rz(-0.60045419) q[3];
sx q[3];
rz(-1.2680935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.101863) q[2];
sx q[2];
rz(-2.1477284) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(-2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(-2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927004) q[0];
sx q[0];
rz(-3.0954439) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(0.13892826) q[1];
sx q[1];
rz(-2.6815092) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28946653) q[0];
sx q[0];
rz(-2.9331101) q[0];
sx q[0];
rz(2.8058488) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87562008) q[2];
sx q[2];
rz(-2.5555829) q[2];
sx q[2];
rz(-0.74967521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13008598) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(1.4873051) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89150724) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(-2.8241983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56269318) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(0.29433027) q[2];
rz(2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-0.99564266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6482553) q[0];
sx q[0];
rz(-2.2264037) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(0.27063453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3753189) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(2.4256698) q[0];
x q[1];
rz(1.0893906) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(-0.35441986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7120413) q[1];
sx q[1];
rz(-1.5263285) q[1];
sx q[1];
rz(-3.0958423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4326914) q[3];
sx q[3];
rz(-2.0376251) q[3];
sx q[3];
rz(2.4469568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3140807) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(-2.6861526) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(-0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6311326) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-0.49490067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1586944) q[0];
sx q[0];
rz(-1.532908) q[0];
sx q[0];
rz(2.057103) q[0];
rz(0.020936326) q[2];
sx q[2];
rz(-0.32792056) q[2];
sx q[2];
rz(2.7431969) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7086664) q[1];
sx q[1];
rz(-2.4927944) q[1];
sx q[1];
rz(1.385958) q[1];
x q[2];
rz(-1.4570191) q[3];
sx q[3];
rz(-0.69996951) q[3];
sx q[3];
rz(-1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0080002) q[2];
sx q[2];
rz(-1.091489) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(1.0333992) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(0.83256759) q[1];
sx q[1];
rz(-1.7201798) q[1];
sx q[1];
rz(-1.7725772) q[1];
rz(-1.9033296) q[2];
sx q[2];
rz(-0.56849545) q[2];
sx q[2];
rz(-0.93760437) q[2];
rz(1.3151863) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
