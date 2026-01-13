function [p_Death]=Death_Probability()
% https://d1wqtxts1xzle7.cloudfront.net/121637207/189-Supplement_1-S69-libre.pdf?1741021642=&response-content-disposition=inline%3B+filename%3DAcute_Measles_Mortality_in_the_United_St.pdf&Expires=1768244579&Signature=CE6FBH9S4ZgvULvLG7QMyyaqtfzLNlVS-iQDr65DBShblUMzd07bg7yNZKuwQDAjppEXrN4z54XcsLss8DHs1XjlySS8CycrbKRqbBFCKwgLFWXgRxgtRlFf7ozSXSBC9esLuE66LAmvvMsakZY3EuteZrpc~Updll~MHoyuD6qLg-DTeMpSZkO5r8Ali6jnupBXl5WWklZN~Rq7NIpwIi6kwZQrhrm7AkZbvL8ZNQcm5kTlukRdNYoAzZb16fSN8N7v82hCbxpU0Z~Yooa1W3idLAQtAWWBGhkehS6K2JwrdrRymZM8Jc8EV1YtALXk3ZfV-6mMCkWfyecSm0CwNQ__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA
% Table 3 NCHS underlying-cause mortality
    p_Death=zeros(1,18);

    p_Death(1)=(25+45)/(9480+18467);
    p_Death(2:4)=27/24114;
    p_Death(5:18)=36/12540;
end