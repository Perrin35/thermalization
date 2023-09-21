OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(-0.54418286) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13574164) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(2.6430921) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40989032) q[2];
sx q[2];
rz(-0.11046834) q[2];
sx q[2];
rz(-0.2875178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4197598) q[1];
sx q[1];
rz(-0.83336035) q[1];
sx q[1];
rz(-1.014773) q[1];
rz(0.55239001) q[3];
sx q[3];
rz(-3.1079709) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.3624181) q[2];
rz(0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966999) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(-1.9702966) q[0];
rz(-0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(1.404095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8453168) q[0];
sx q[0];
rz(-0.8527841) q[0];
sx q[0];
rz(2.5677471) q[0];
rz(-pi) q[1];
rz(-0.48268433) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(2.5978616) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5215056) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(0.28316811) q[1];
rz(-1.4349798) q[3];
sx q[3];
rz(-2.0471626) q[3];
sx q[3];
rz(-1.7686896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(0.94397604) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.336162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29207644) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-2.3220093) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(-0.56366411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.992313) q[0];
sx q[0];
rz(-0.57846071) q[0];
sx q[0];
rz(3.0119386) q[0];
x q[1];
rz(-1.7965505) q[2];
sx q[2];
rz(-2.1588237) q[2];
sx q[2];
rz(-0.079344582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.015805294) q[1];
sx q[1];
rz(-1.8965221) q[1];
sx q[1];
rz(-0.9822407) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1590957) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(1.8713649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(1.0602661) q[2];
rz(3.0660196) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(-2.9484205) q[0];
rz(1.5974143) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(2.4838122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5415216) q[0];
sx q[0];
rz(-2.9642448) q[0];
sx q[0];
rz(0.17139165) q[0];
x q[1];
rz(2.8305956) q[2];
sx q[2];
rz(-1.9823091) q[2];
sx q[2];
rz(2.3222773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12280497) q[1];
sx q[1];
rz(-1.0446761) q[1];
sx q[1];
rz(-0.29463525) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85122078) q[3];
sx q[3];
rz(-0.38098601) q[3];
sx q[3];
rz(-1.756543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(-0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(-0.62686282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1889362) q[0];
sx q[0];
rz(-3.1075826) q[0];
sx q[0];
rz(-1.0556428) q[0];
rz(-pi) q[1];
rz(-2.3847694) q[2];
sx q[2];
rz(-2.7177817) q[2];
sx q[2];
rz(-0.28048453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53351952) q[1];
sx q[1];
rz(-2.7058209) q[1];
sx q[1];
rz(2.3418952) q[1];
x q[2];
rz(0.12368006) q[3];
sx q[3];
rz(-1.9607753) q[3];
sx q[3];
rz(-2.9434413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(1.0181001) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-2.0904082) q[1];
sx q[1];
rz(-0.32454023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6231461) q[0];
sx q[0];
rz(-2.0707957) q[0];
sx q[0];
rz(-2.6536233) q[0];
x q[1];
rz(2.2424477) q[2];
sx q[2];
rz(-2.2567344) q[2];
sx q[2];
rz(-1.9859973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2734087) q[1];
sx q[1];
rz(-1.5955828) q[1];
sx q[1];
rz(2.1048057) q[1];
rz(-1.8683897) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(-0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4679608) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(-2.0868789) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(1.2021525) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5613272) q[0];
sx q[0];
rz(-2.1939477) q[0];
sx q[0];
rz(-2.6448963) q[0];
x q[1];
rz(0.41263327) q[2];
sx q[2];
rz(-0.79421439) q[2];
sx q[2];
rz(0.40976322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.051604328) q[1];
sx q[1];
rz(-3.0172536) q[1];
sx q[1];
rz(-2.405637) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2644516) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(0.84623519) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071035944) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.6092469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198062) q[0];
sx q[0];
rz(-2.0397423) q[0];
sx q[0];
rz(1.8437587) q[0];
rz(-0.74560994) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(-1.2127884) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6541877) q[1];
sx q[1];
rz(-2.020917) q[1];
sx q[1];
rz(-0.022189157) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5395457) q[3];
sx q[3];
rz(-1.6685408) q[3];
sx q[3];
rz(2.5412113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9939076) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.9130075) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(-1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.5478741) q[0];
sx q[0];
rz(-2.6429208) q[0];
sx q[0];
rz(-0.6643995) q[0];
rz(1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(1.857035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212095) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(-2.3032805) q[0];
x q[1];
rz(-2.8753488) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(-3.0014696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38332332) q[1];
sx q[1];
rz(-2.3136534) q[1];
sx q[1];
rz(2.5187056) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5981538) q[3];
sx q[3];
rz(-1.1809071) q[3];
sx q[3];
rz(2.011812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(0.55076304) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-3.0974467) q[0];
rz(-1.6126397) q[1];
sx q[1];
rz(-1.2395369) q[1];
sx q[1];
rz(-2.3619161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6990307) q[0];
sx q[0];
rz(-2.583722) q[0];
sx q[0];
rz(2.872422) q[0];
x q[1];
rz(1.0829955) q[2];
sx q[2];
rz(-2.3239845) q[2];
sx q[2];
rz(1.4873193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1624987) q[1];
sx q[1];
rz(-1.7531803) q[1];
sx q[1];
rz(2.6136293) q[1];
rz(-pi) q[2];
rz(0.78124222) q[3];
sx q[3];
rz(-1.2282073) q[3];
sx q[3];
rz(2.4398746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4298532) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(1.7977057) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-2.8316108) q[2];
sx q[2];
rz(-1.1520755) q[2];
sx q[2];
rz(0.68785695) q[2];
rz(2.963831) q[3];
sx q[3];
rz(-2.4958785) q[3];
sx q[3];
rz(2.6444825) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];