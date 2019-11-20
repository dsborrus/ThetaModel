Hello Greg,

I've put everything you'll need into one folder to try and simplify everything.

The code to run, and have it just GO, is main.m

main.m assigns all the parameters, and script options, such as which plots to make, whether to do ablations or which type of network to use etc. etc. ...

But then, main.m calls the real heavy weight script (function) simulate_v1.m. simulate_v1.m is where the simulation actually takes place. You shouldn't have to touch it much, except maybe for plotting. Personally, I actually have simulate_v1.m CALL another script to do plotting just to keep it seperate.

I've set it up so main.m should should just run no problem. You'll have to play around in simulate_v1.m in order to see the names of the state variables you're interested in.

Let me know if you need any help!

-DB 190715


