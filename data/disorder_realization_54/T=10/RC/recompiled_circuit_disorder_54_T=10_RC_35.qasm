OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(3.830885) q[0];
sx q[0];
rz(9.7552714) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56624428) q[0];
sx q[0];
rz(-0.69501221) q[0];
sx q[0];
rz(1.1027176) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.032151392) q[2];
sx q[2];
rz(-1.8841779) q[2];
sx q[2];
rz(-1.9762447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6735437) q[1];
sx q[1];
rz(-1.1357726) q[1];
sx q[1];
rz(-0.97462868) q[1];
rz(0.96592824) q[3];
sx q[3];
rz(-2.5336207) q[3];
sx q[3];
rz(-2.0917497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(0.7888166) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(-0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5641862) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(-2.2400411) q[0];
x q[1];
rz(-3.0099478) q[2];
sx q[2];
rz(-2.3028214) q[2];
sx q[2];
rz(1.6694348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20827046) q[1];
sx q[1];
rz(-1.9083438) q[1];
sx q[1];
rz(-2.7949105) q[1];
rz(-pi) q[2];
x q[2];
rz(1.968859) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(1.845899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(2.3201578) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-2.6170513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6150104) q[0];
sx q[0];
rz(-0.16631642) q[0];
sx q[0];
rz(1.7517356) q[0];
rz(-1.8345941) q[2];
sx q[2];
rz(-2.2321777) q[2];
sx q[2];
rz(1.7922572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4265392) q[1];
sx q[1];
rz(-1.5476777) q[1];
sx q[1];
rz(3.0798562) q[1];
rz(-pi) q[2];
rz(-1.8936162) q[3];
sx q[3];
rz(-0.4399235) q[3];
sx q[3];
rz(0.80252121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.4952205) q[2];
rz(1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33048531) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(0.088949732) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8331497) q[0];
sx q[0];
rz(-1.8652548) q[0];
sx q[0];
rz(-1.740728) q[0];
x q[1];
rz(-1.7961411) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.00011132414) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(-2.6898756) q[1];
rz(-pi) q[2];
rz(2.8702909) q[3];
sx q[3];
rz(-2.6415083) q[3];
sx q[3];
rz(0.015319583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(-2.518667) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(0.46666551) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899807) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.4978283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72901112) q[0];
sx q[0];
rz(-0.81775613) q[0];
sx q[0];
rz(1.9304995) q[0];
rz(0.79767144) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(-2.5728512) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6240546) q[1];
sx q[1];
rz(-2.2741286) q[1];
sx q[1];
rz(-1.5348977) q[1];
x q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-1.0609396) q[3];
sx q[3];
rz(1.7061403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(1.5779457) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-2.5031228) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88046056) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(2.3418505) q[0];
x q[1];
rz(-3.0827423) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(-2.5318052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1286436) q[1];
sx q[1];
rz(-1.1638906) q[1];
sx q[1];
rz(-3.0055771) q[1];
rz(-pi) q[2];
rz(1.6834429) q[3];
sx q[3];
rz(-0.23922353) q[3];
sx q[3];
rz(0.84157543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0662213) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(-0.071468778) q[2];
rz(1.4525157) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(-0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(-1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(-1.0095899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1560695) q[0];
sx q[0];
rz(-2.0445604) q[0];
sx q[0];
rz(1.9457293) q[0];
x q[1];
rz(1.3426443) q[2];
sx q[2];
rz(-1.9462799) q[2];
sx q[2];
rz(2.0827039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13751444) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(0.60933463) q[1];
rz(0.76241242) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(-1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1040161) q[0];
sx q[0];
rz(-0.74066478) q[0];
sx q[0];
rz(-0.91233493) q[0];
rz(-pi) q[1];
rz(-0.93034805) q[2];
sx q[2];
rz(-1.010251) q[2];
sx q[2];
rz(0.7330187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6726482) q[1];
sx q[1];
rz(-0.69744195) q[1];
sx q[1];
rz(-1.3717321) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29856155) q[3];
sx q[3];
rz(-2.7136554) q[3];
sx q[3];
rz(-0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(1.0347962) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(2.8544193) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(0.33219355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848541) q[0];
sx q[0];
rz(-0.85002725) q[0];
sx q[0];
rz(2.7194354) q[0];
x q[1];
rz(-0.046499649) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(0.38682129) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13739535) q[1];
sx q[1];
rz(-1.9837556) q[1];
sx q[1];
rz(1.2910299) q[1];
x q[2];
rz(-2.6094749) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(-1.2906637) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(0.97958952) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-2.8870781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33728889) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(-1.8810349) q[0];
rz(0.48760957) q[2];
sx q[2];
rz(-1.3541823) q[2];
sx q[2];
rz(0.85607869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0490129) q[1];
sx q[1];
rz(-1.4900467) q[1];
sx q[1];
rz(-2.8156274) q[1];
rz(-2.7032095) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(2.5509978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(1.7808328) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(-2.5361983) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1142674) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(2.7813773) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-2.0965626) q[2];
sx q[2];
rz(-1.7797995) q[2];
sx q[2];
rz(-1.9893653) q[2];
rz(1.249282) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
