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
rz(-1.2000097) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8631247) q[0];
sx q[0];
rz(-1.5764563) q[0];
sx q[0];
rz(1.6186884) q[0];
rz(2.0179022) q[2];
sx q[2];
rz(-1.8748218) q[2];
sx q[2];
rz(2.6927039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0389509) q[1];
sx q[1];
rz(-2.0248374) q[1];
sx q[1];
rz(-2.7137043) q[1];
rz(-pi) q[2];
rz(1.5267738) q[3];
sx q[3];
rz(-1.9485954) q[3];
sx q[3];
rz(1.6954741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.6480185) q[2];
sx q[2];
rz(2.8519894) q[2];
rz(-2.2662207) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-0.059710596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5263379) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-2.7976024) q[0];
rz(-0.084331766) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(-1.7864236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8665168) q[0];
sx q[0];
rz(-1.6356902) q[0];
sx q[0];
rz(-2.1823723) q[0];
rz(1.6927035) q[2];
sx q[2];
rz(-1.4093471) q[2];
sx q[2];
rz(0.09300692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.018108798) q[1];
sx q[1];
rz(-1.5068441) q[1];
sx q[1];
rz(2.0204087) q[1];
rz(-pi) q[2];
rz(-1.1739028) q[3];
sx q[3];
rz(-0.19860425) q[3];
sx q[3];
rz(1.7484776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8460059) q[2];
sx q[2];
rz(-0.83335525) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(2.6387571) q[3];
sx q[3];
rz(-0.42481315) q[3];
sx q[3];
rz(-2.1311549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-0.53776598) q[0];
sx q[0];
rz(-0.64087254) q[0];
rz(-2.3928941) q[1];
sx q[1];
rz(-2.1332032) q[1];
sx q[1];
rz(-1.0650939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051620313) q[0];
sx q[0];
rz(-1.9301156) q[0];
sx q[0];
rz(-3.1348455) q[0];
rz(-pi) q[1];
rz(-2.2674019) q[2];
sx q[2];
rz(-1.0993996) q[2];
sx q[2];
rz(-0.092560571) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47485891) q[1];
sx q[1];
rz(-1.534728) q[1];
sx q[1];
rz(-0.95400793) q[1];
rz(-1.8530811) q[3];
sx q[3];
rz(-2.3989587) q[3];
sx q[3];
rz(1.7395072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(2.4242145) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-2.513956) q[3];
sx q[3];
rz(-0.29754105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3142969) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(-2.1266134) q[1];
sx q[1];
rz(-0.29622886) q[1];
sx q[1];
rz(3.1304741) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7751559) q[0];
sx q[0];
rz(-1.9754793) q[0];
sx q[0];
rz(1.0767656) q[0];
x q[1];
rz(-2.501802) q[2];
sx q[2];
rz(-0.94546972) q[2];
sx q[2];
rz(1.7196136) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0032469) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(-0.37133118) q[1];
rz(-pi) q[2];
rz(1.3783781) q[3];
sx q[3];
rz(-1.8835861) q[3];
sx q[3];
rz(-2.6421412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6461688) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(-2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-0.88808131) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11113142) q[0];
sx q[0];
rz(-0.24065329) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(-1.8146851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3834284) q[0];
sx q[0];
rz(-0.82793068) q[0];
sx q[0];
rz(-1.3616256) q[0];
x q[1];
rz(2.5577776) q[2];
sx q[2];
rz(-2.2540255) q[2];
sx q[2];
rz(0.77606397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8718308) q[1];
sx q[1];
rz(-2.7360536) q[1];
sx q[1];
rz(0.62015066) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29019659) q[3];
sx q[3];
rz(-1.962933) q[3];
sx q[3];
rz(-2.7249667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1422687) q[2];
sx q[2];
rz(-2.6596255) q[2];
sx q[2];
rz(1.8959321) q[2];
rz(1.2549531) q[3];
sx q[3];
rz(-0.84147036) q[3];
sx q[3];
rz(-2.3033223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44928837) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(2.4601049) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(2.025827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9003446) q[0];
sx q[0];
rz(-1.1660518) q[0];
sx q[0];
rz(2.3099398) q[0];
rz(-pi) q[1];
rz(0.90066465) q[2];
sx q[2];
rz(-1.2256983) q[2];
sx q[2];
rz(-2.5010441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15397729) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(-1.5495367) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9168105) q[3];
sx q[3];
rz(-2.4596679) q[3];
sx q[3];
rz(3.0608321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21268022) q[2];
sx q[2];
rz(-1.8447515) q[2];
sx q[2];
rz(0.35432717) q[2];
rz(2.8220693) q[3];
sx q[3];
rz(-1.1525681) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(-1.1122423) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(2.5792714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4903957) q[0];
sx q[0];
rz(-2.2389452) q[0];
sx q[0];
rz(-3.1227123) q[0];
rz(-2.8096335) q[2];
sx q[2];
rz(-1.6662285) q[2];
sx q[2];
rz(1.7711668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9971804) q[1];
sx q[1];
rz(-1.4257396) q[1];
sx q[1];
rz(-0.1952862) q[1];
x q[2];
rz(0.9741707) q[3];
sx q[3];
rz(-1.6430292) q[3];
sx q[3];
rz(2.7330287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.101863) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(2.9240821) q[3];
sx q[3];
rz(-2.1931931) q[3];
sx q[3];
rz(2.8229009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927004) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(0.39644077) q[0];
rz(-0.13892826) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28946653) q[0];
sx q[0];
rz(-0.20848256) q[0];
sx q[0];
rz(0.33574386) q[0];
rz(-pi) q[1];
rz(1.0993768) q[2];
sx q[2];
rz(-1.9328914) q[2];
sx q[2];
rz(0.21381703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8926881) q[1];
sx q[1];
rz(-0.78033328) q[1];
sx q[1];
rz(-3.0569539) q[1];
rz(0.89150724) q[3];
sx q[3];
rz(-1.659698) q[3];
sx q[3];
rz(0.31739435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5788995) q[2];
sx q[2];
rz(-1.0572628) q[2];
sx q[2];
rz(2.8472624) q[2];
rz(-1.1307905) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(-0.99564266) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482553) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-2.1954779) q[1];
sx q[1];
rz(-0.39603907) q[1];
sx q[1];
rz(2.8709581) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3753189) q[0];
sx q[0];
rz(-0.047485504) q[0];
sx q[0];
rz(0.71592285) q[0];
rz(-2.052202) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(-0.35441986) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86079019) q[1];
sx q[1];
rz(-1.5250912) q[1];
sx q[1];
rz(1.526282) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6587894) q[3];
sx q[3];
rz(-0.82595347) q[3];
sx q[3];
rz(0.39289075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82751194) q[2];
sx q[2];
rz(-0.22958799) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(2.3296302) q[3];
sx q[3];
rz(-0.88192314) q[3];
sx q[3];
rz(0.5493831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(-0.73927885) q[0];
rz(-0.23070681) q[1];
sx q[1];
rz(-2.6735327) q[1];
sx q[1];
rz(-2.646692) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65942818) q[0];
sx q[0];
rz(-0.48766252) q[0];
sx q[0];
rz(-1.6517261) q[0];
rz(-pi) q[1];
rz(0.020936326) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(0.39839572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.478185) q[1];
sx q[1];
rz(-2.2067398) q[1];
sx q[1];
rz(-3.0031167) q[1];
rz(-pi) q[2];
rz(-0.87402113) q[3];
sx q[3];
rz(-1.6439983) q[3];
sx q[3];
rz(-0.30833581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13359244) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(2.8137394) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-2.8938507) q[3];
sx q[3];
rz(1.0333992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1223758) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(-0.83256759) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-2.1140425) q[2];
sx q[2];
rz(-1.3941358) q[2];
sx q[2];
rz(0.91640581) q[2];
rz(-1.8264063) q[3];
sx q[3];
rz(-1.8288463) q[3];
sx q[3];
rz(2.0127206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];